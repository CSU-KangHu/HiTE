import os
import yaml
import logging
import argparse
import numpy as np
import pandas as pd
from tqdm import tqdm
import torch
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F
from torch.nn.parallel import DataParallel
from torch.utils.data import TensorDataset, DataLoader, random_split
from sklearn.utils import shuffle
from sklearn.metrics import classification_report, precision_score, accuracy_score, recall_score, f1_score, roc_auc_score, roc_curve, auc, precision_recall_curve, average_precision_score
from model import CNNCAT
from sklearn.calibration import CalibratedClassifierCV
from torch.optim.lr_scheduler import OneCycleLR, CosineAnnealingWarmRestarts
import joblib
from sklearn.linear_model import LogisticRegression
from focal_loss import FocalLoss
from sklearn.model_selection import KFold, StratifiedKFold

# Set a global seed for reproducibility
SEED = 2024
np.random.seed(SEED)
torch.manual_seed(SEED)
if torch.cuda.is_available():
    torch.cuda.manual_seed_all(SEED)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False

class LoggerMaker:
    """A helper class to create and configure a logger."""
    def __init__(self, out_dir, log_name='global', log_file=None, level=logging.DEBUG):
        self.out_dir = out_dir
        self.logger = None
        self.log_file = log_file
        self.log_name = log_name
        self.level = level
        self._make_logger()
    
    def _make_logger(self):
        """Initializes the logger with specified handlers and formatters."""
        logger = logging.getLogger(self.log_name)
        logger.setLevel(self.level)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

        # Console handler
        ch = logging.StreamHandler()
        ch.setLevel(self.level)
        ch.setFormatter(formatter)
        logger.addHandler(ch)

        # File handler
        if self.log_file:
            fh = logging.FileHandler(os.path.join(self.out_dir, self.log_file))
            fh.setLevel(self.level)
            fh.setFormatter(formatter)
            logger.addHandler(fh)

        self.logger = logger

class CNNTrainer:
    """A class to handle the training, validation, and testing of the CNN model."""
    def __init__(self, num_epochs: int, learning_rate: float, l2_weight: float, batch_size: int, device_id, logger):
        self.output_dir = './'
        self.num_epochs = num_epochs
        self.learning_rate = learning_rate
        self.l2_weight = l2_weight
        self.batch_size = batch_size
        self.device_id = device_id
        # Set the device for training (GPU or CPU)
        self.device = torch.device(f'cuda:{self.device_id[0]}' if torch.cuda.is_available() else 'cpu')
        # Set up logger
        self.logger = logger
      
        self.logits_list = []
        self.labels_list = []

        self.train_loss = []
        self.val_loss = []

        self.logger.info("CNNTrainer initialized successfully")

    def _save_checkpoint(self, model: nn.Module, optimizer: optim.Optimizer, model_name, 
                        epoch: int, metrics: dict, is_best: bool = False, fold: int = None):
        """
        Saves a model checkpoint.
        
        Args:
            model: The model to save.
            optimizer: The optimizer state.
            epoch: The current epoch number.
            metrics: A dictionary of performance metrics.
            is_best: A flag to indicate if this is the best model so far.
            fold: The current cross-validation fold number.
        """
        checkpoint = {
            'epoch': epoch,
            'model_state_dict': model.state_dict(),
            'optimizer_state_dict': optimizer.state_dict(),
            'metrics': metrics,
            'fold': fold
        }
        
        # Save the latest checkpoint
        checkpoint_path = os.path.join(
            self.output_dir,
            model_name + '.pth'
        )
        torch.save(checkpoint, checkpoint_path)
        self.logger.info(f"Saved checkpoint at {checkpoint_path}")

    def _setup_scheduler(self, optimizer, num_epochs, train_loader):
        """Sets up the learning rate scheduler."""
        scheduler_type = 'cosine'
        
        if scheduler_type == 'onecycle':
            scheduler = OneCycleLR(
                optimizer,
                max_lr=0.01,
                epochs=num_epochs,
                steps_per_epoch=len(train_loader),
                pct_start=0.3,
                anneal_strategy='cos'
            )
        elif scheduler_type == 'cosine':
            scheduler = CosineAnnealingWarmRestarts(
                optimizer,
                T_0=5,
                T_mult=2,
                eta_min=1e-6
            )
        else:
            raise ValueError(f"Unknown scheduler type: {scheduler_type}")
            
        return scheduler

    def _compute_metrics(self, labels: np.ndarray, probs: np.ndarray, predictions: np.ndarray) -> dict:
        """
        Computes classification metrics.
        
        Args:
            labels: True labels.
            probs: Predicted probabilities for the positive class.
            predictions: Predicted class labels.
            
        Returns:
            dict: A dictionary of classification metrics.
        """
        # Binarize predictions based on a 0.5 threshold
        predictions = (probs > 0.5).astype(int)
        accuracy = accuracy_score(labels, predictions)
        precision = precision_score(
            labels, predictions, average='macro', zero_division=0
        )
        recall = recall_score(
            labels, predictions, average='macro', zero_division=0
        )
        f1 = f1_score(labels, predictions, average='macro', zero_division=0)

        roc_auc = roc_auc_score(labels, probs)
        pr_auc = average_precision_score(labels, probs)
        
        return {
            'accuracy': accuracy,
            'precision': precision,
            'recall': recall,
            'f1': f1,
            'roc_auc': roc_auc,
            'pr_auc': pr_auc
        }

    def train_epoch(self, model: nn.Module, optimizer: optim.Optimizer, scheduler: optim.lr_scheduler._LRScheduler,
                   criterion: nn.Module, train_loader: DataLoader,
                   epoch: int) -> tuple:
        """
        Trains the model for one epoch.
        
        Args:
            model: The neural network model.
            optimizer: The optimizer.
            scheduler: The learning rate scheduler.
            criterion: The loss function.
            train_loader: The data loader for training data.
            epoch: The current epoch number.
            
        Returns:
            A tuple containing average loss and a dictionary of metrics.
        """
        model.train()
        total_loss = 0
        epoch_probs = []
        epoch_labels = []
        
        for batch_x, batch_x_kmer, batch_y in train_loader:
            try:
                # Prepare data and move to device
                batch_x = batch_x.to(self.device)
                batch_x_kmer = batch_x_kmer.to(self.device)
                batch_y = batch_y.to(self.device)
                
                # Forward pass
                optimizer.zero_grad()
                logits, l2_reg = model(batch_x, batch_x_kmer)
                loss = criterion(logits, batch_y)
                
                # Backward pass and optimization
                loss.backward()
                optimizer.step()
                scheduler.step()
                
                # Record predictions and labels for metric calculation
                probs = F.softmax(logits, dim=1)[:, 1].cpu().detach().numpy()
                epoch_probs.extend(probs)
                epoch_labels.extend(batch_y.cpu().numpy())
                total_loss += loss.item()
                
            except Exception as e:
                self.logger.error(f"Error in training batch: {str(e)}")
                raise
                
        # Compute metrics for the epoch
        avg_loss = total_loss / len(train_loader)
        self.train_loss.append(avg_loss)
        metrics = self._compute_metrics(
            np.array(epoch_labels), np.array(epoch_probs), None
        )
       
        return avg_loss, metrics

    def validate(self, model: nn.Module, criterion: nn.Module,
                val_loader: DataLoader) -> tuple:
        """
        Evaluates the model on the validation set.
        
        Args:
            model: The neural network model.
            criterion: The loss function.
            val_loader: The data loader for validation data.
            
        Returns:
            A tuple containing average loss and a dictionary of metrics.
        """
        model.eval()
        total_loss = 0
        epoch_probs = []
        epoch_labels = []
        
        with torch.no_grad():
            for batch_x, batch_x_kmer, batch_y in val_loader:
                try:
                    batch_x = batch_x.to(self.device)
                    batch_x_kmer = batch_x_kmer.to(self.device)
                    batch_y = batch_y.to(self.device)
                    
                    logits, l2_reg = model(batch_x, batch_x_kmer)
                    loss = criterion(logits, batch_y)
                    total_loss += loss.item() 
                    
                    probs = F.softmax(logits, dim=1)[:, 1].cpu().detach().numpy()
                    epoch_probs.extend(probs)
                    epoch_labels.append(batch_y.cpu().numpy())
                    
                except Exception as e:
                    self.logger.error(f"Error in validation batch: {str(e)}")
                    raise
        
        # Aggregate results from all batches
        prob_val = np.concatenate(epoch_probs)
        y_val = np.concatenate(epoch_labels)
        
        # Compute metrics
        avg_loss = total_loss / len(val_loader)
        self.val_loss.append(avg_loss)
        metrics = self._compute_metrics(y_val, prob_val, None)
        
        return avg_loss, metrics

    def test(self, model: nn.Module, criterion: nn.Module,
             test_loader: DataLoader) -> tuple:
        """
        Evaluates the model on the test set.
        
        Args:
            model: The neural network model.
            criterion: The loss function.
            test_loader: The data loader for test data.
            
        Returns:
            A tuple containing average loss and a dictionary of metrics.
        """
        model.eval()
        total_loss = 0
        epoch_probs = []
        epoch_labels = []
        
        with torch.no_grad():
            for batch_x, batch_x_kmer, batch_y in tqdm(test_loader, desc='Testing'):
                try:
                    batch_x = batch_x.to(self.device)
                    batch_x_kmer = batch_x_kmer.to(self.device)
                    batch_y = batch_y.to(self.device)
                    
                    logits, l2_reg = model(batch_x, batch_x_kmer)
                    loss = criterion(logits, batch_y)
                    total_loss += loss.item()
                    
                    probs = F.softmax(logits, dim=1)[:, 1].cpu().detach().numpy()
                    epoch_probs.extend(probs)
                    epoch_labels.append(batch_y.cpu().numpy())
                                     
                except Exception as e:
                    self.logger.error(f"Error in test batch: {str(e)}")
                    raise
        
        # Aggregate results from all batches
        probs_val = np.concatenate(epoch_probs)
        y_val = np.concatenate(epoch_labels)
       
        # Compute metrics
        avg_loss = total_loss / len(test_loader)
        metrics = self._compute_metrics(y_val, probs_val, None)

        return avg_loss, metrics

    def train(self, train_dataset, test_dataset):
        """Executes the full training and evaluation pipeline using cross-validation."""
        try:
            # Get training parameters
            num_epoch = self.num_epochs
            lr = self.learning_rate
            l2 = self.l2_weight
            bs = self.batch_size
            n_folds = 5
            
            all_metrics = []
            random_seeds = [2024, 2025, 2026, 2027, 2028]
            labels = train_dataset.tensors[2]

            # Repeat 5-fold cross-validation 5 times with different random seeds
            for i in range(5):
                self.logger.info(f"\n CV Run: {i+1}")
                fold_metrics = []
                # Create a stratified K-fold splitter
                skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=random_seeds[i])
                
                for fold, (train_idx, val_idx) in enumerate(skf.split(torch.zeros(len(labels)), labels)):
                    self.train_loss = []
                    self.val_loss = []
                    self.best_loss = float('inf')
                    self.logger.info(f"\nTraining Fold {fold + 1}/{n_folds}")

                    # Create data subsets for training and validation
                    train_subset = torch.utils.data.Subset(train_dataset, train_idx)
                    val_subset = torch.utils.data.Subset(train_dataset, val_idx)

                    # Create data loaders
                    train_loader = DataLoader(train_subset, batch_size=bs, shuffle=True, drop_last=True)
                    val_loader = DataLoader(val_subset, batch_size=bs)

                    # Initialize model
                    model = CNNCAT().to(self.device)

                    # Use DataParallel for multi-GPU training
                    if len(self.device_id) > 1 and torch.cuda.device_count() > 1:
                        self.device_id = [int(dev) for dev in self.device_id]
                        self.logger.info(f"Using GPUs: {self.device_id}")
                        model = DataParallel(model, device_ids=self.device_id)
                        
                    # Set loss function
                    criterion = FocalLoss()

                    # Set optimizer and scheduler
                    optimizer = optim.AdamW(
                        model.parameters(),
                        lr=lr,
                        weight_decay=l2
                    )
                    scheduler = self._setup_scheduler(optimizer, num_epoch, train_loader)

                    # Track the best metrics for the current fold
                    fold_best_metrics = {'pr_auc': 0}

                    # Training loop
                    for epoch in range(num_epoch):
                        train_loss, train_metrics = self.train_epoch(
                            model, optimizer, scheduler, criterion, train_loader, epoch
                        )
                        val_loss, val_metrics = self.validate(
                            model, criterion, val_loader
                        )
                        if (epoch + 1) % 5 == 0:
                            self.logger.info(
                                f"\nFold {fold + 1}, Epoch {epoch+1}/{num_epoch}, Train Loss: {train_loss:.5f}, "
                                f"Val Loss: {val_loss:.5f}, Val PR AUC: {val_metrics['pr_auc']:.5f}, Val ROC AUC: {val_metrics['roc_auc']:.5f}"
                            )
                    
                        # Update best metrics for the fold based on validation PR AUC
                        if val_metrics['pr_auc'] > fold_best_metrics['pr_auc']:
                            fold_best_metrics.update(val_metrics)

                    # Evaluate the final model of the fold on the test set
                    test_dataloader = DataLoader(test_dataset, batch_size=bs)
                    avg_loss, test_metrics = self.test(model, criterion, test_dataloader)
                    self.logger.info(
                        f"\nFold {fold + 1} Test Results -> Loss: {avg_loss:.5f}, PR AUC: {test_metrics['pr_auc']:.5f}, ROC AUC: {test_metrics['roc_auc']:.5f}\n"
                    )
                    fold_metrics.append(test_metrics)

                    # Save a checkpoint of the final model for this fold
                    self._save_checkpoint(
                        model, optimizer, f'checkpoint_run{i+1}_fold{fold+1}', epoch, fold_best_metrics, fold=fold
                    )

                # Calculate and log the average metrics across all folds for the current run
                avg_metrics = {
                    metric: np.mean([fold[metric] for fold in fold_metrics])
                    for metric in ['accuracy', 'precision', 'recall', 'f1', 'roc_auc', 'pr_auc']
                }
                all_metrics.append(avg_metrics)
                self.logger.info(
                    f"\nCV Run {i+1} Average Results ({n_folds} folds):\n"
                    f"Accuracy: {avg_metrics['accuracy']:.5f} ± {np.std([f['accuracy'] for f in fold_metrics]):.5f}\n"
                    f"Precision: {avg_metrics['precision']:.5f} ± {np.std([f['precision'] for f in fold_metrics]):.5f}\n"
                    f"Recall: {avg_metrics['recall']:.5f} ± {np.std([f['recall'] for f in fold_metrics]):.5f}\n"
                    f"F1: {avg_metrics['f1']:.5f} ± {np.std([f['f1'] for f in fold_metrics]):.5f}\n"
                    f"ROC AUC: {avg_metrics['roc_auc']:.5f} ± {np.std([f['roc_auc'] for f in fold_metrics]):.5f}\n"
                    f"PR AUC: {avg_metrics['pr_auc']:.5f} ± {np.std([f['pr_auc'] for f in fold_metrics]):.5f}\n"
                )

            return all_metrics  

        except Exception as e:
            self.logger.error(f"Training failed: {str(e)}")
            raise

def load_full_dataset(input_dir, logger):
    """
    Loads and preprocesses the full dataset.
    
    Returns:
        A tuple of TensorDataset for training and testing.
    """
    def _normalize_features(features: torch.Tensor) -> torch.Tensor:
        """Normalizes features using L2 normalization."""
        try:
            return F.normalize(features.float(), p=2, dim=1)
        except Exception as e:
            print(f"Feature normalization failed: {str(e)}")
            raise

    try:
        # Load training data
        train_img_features = _normalize_features(torch.load(os.path.join(input_dir, 'train_img_features.pt')))
        train_freq_features = _normalize_features(torch.load(os.path.join(input_dir, 'train_freq_features.pt')))
        train_labels = torch.load(os.path.join(input_dir, 'train_labels.pt'))
        
        # Load testing data
        test_img_features = _normalize_features(torch.load(os.path.join(input_dir, 'test_img_features.pt')))
        test_freq_features = _normalize_features(torch.load(os.path.join(input_dir, 'test_freq_features.pt')))
        test_labels = torch.load(os.path.join(input_dir, 'test_labels.pt')).long()
        
        logger.info(f'Train label distribution:\n{pd.Series(train_labels.numpy()).value_counts()}')
        logger.info(f'Test label distribution:\n{pd.Series(test_labels.numpy()).value_counts()}')

        # Create datasets
        train_dataset = TensorDataset(train_img_features, train_freq_features, train_labels)
        test_dataset = TensorDataset(test_img_features, test_freq_features, test_labels)
                
        return train_dataset, test_dataset
        
    except Exception as e:
        print(f"Full dataset loading failed: {str(e)}")
        raise

def train_best(train_dataset, test_dataset, device_id, logger):
    """Initializes the trainer and starts the training process with the best hyperparameters."""
    trainer = CNNTrainer(
        num_epochs=14,
        learning_rate=0.001,
        l2_weight=0.1,
        batch_size=64,
        device_id=device_id,
        logger=logger
    )

    all_metrics = trainer.train(train_dataset, test_dataset)

    # Log final results across all 5 runs
    final_avg_metrics = {
        metric: np.mean([run[metric] for run in all_metrics])
        for metric in ['accuracy', 'precision', 'recall', 'f1', 'roc_auc', 'pr_auc']
    }

    logger.info(
        f"\n--- Overall Average Results (5 runs of 5-fold CV) ---\n"
        f"Accuracy: {final_avg_metrics['accuracy']:.5f} ± {np.std([f['accuracy'] for f in all_metrics]):.5f}\n"
        f"Precision: {final_avg_metrics['precision']:.5f} ± {np.std([f['precision'] for f in all_metrics]):.5f}\n"
        f"Recall: {final_avg_metrics['recall']:.5f} ± {np.std([f['recall'] for f in all_metrics]):.5f}\n"
        f"F1: {final_avg_metrics['f1']:.5f} ± {np.std([f['f1'] for f in all_metrics]):.5f}\n"
        f"ROC AUC: {final_avg_metrics['roc_auc']:.5f} ± {np.std([f['roc_auc'] for f in all_metrics]):.5f}\n"
        f"PR AUC: {final_avg_metrics['pr_auc']:.5f} ± {np.std([f['pr_auc'] for f in all_metrics]):.5f}\n"
    )
    return final_avg_metrics['pr_auc']

def main():
    input_dir = './' # Your_Feature_Directory_Here
    logger = LoggerMaker('./', 'tune', 'tune.log', logging.INFO).logger
    device_id = [3]
    try:
        train_dataset, test_dataset = load_full_dataset(input_dir, logger)
        train_best(train_dataset, test_dataset, device_id, logger)
    except Exception as e:
        logger.error(f"Training failed: {str(e)}")
        raise

if __name__ == '__main__':
    main()