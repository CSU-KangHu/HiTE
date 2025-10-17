import os
import torch
import logging
import numpy as np
import torch.nn.functional as F
from torch.utils.data import TensorDataset, DataLoader
from typing import List, Tuple, Dict, Union
from tqdm import tqdm
from hybridLTR_model import CNNCAT
import argparse
import pandas as pd

class CNNPredictor:
    def __init__(self, model_path: str, threshold: float, device: str = None, logger=None):
        """
        initialize the predictor
        
        Args:
            model_path: model file path
            device: device type ('cuda' or 'cpu')
        """
        if logger is None:
            self._setup_logging()
        else:
            self.logger = logger

        self.device = device if device else \
                     ('cuda' if torch.cuda.is_available() else 'cpu')
        self.logger.info(f"Using device: {self.device}")
        
        # load model
        self.model = self._load_model(model_path)
        self.model.eval()
        self.threshold = threshold
        print('threshold:', self.threshold)

    def _setup_logging(self):
        """set up logging system"""
        self.logger = logging.getLogger('HybirdLTR')
        self.logger.setLevel(logging.INFO)
        
        # add console handler
        ch = logging.StreamHandler()
        ch.setLevel(logging.DEBUG)
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )
        ch.setFormatter(formatter)
        self.logger.addHandler(ch)

    def _load_model(self, model_path: str) -> torch.nn.Module:
        """
        load the trained model
        
        Args:
            model_path: model file path
            
        Returns:
            torch.nn.Module: loaded model
        """
        try:
            self.logger.info(f"Loading model from {model_path}")
            checkpoint = torch.load(model_path, map_location=self.device)
            
            # initialize the model
            model = CNNCAT().to(self.device)
            
            # if the model is saved with DataParallel, need to handle it
            if list(checkpoint['model_state_dict'].keys())[0].startswith('module.'):
                if self.device == 'cuda':
                    model = torch.nn.DataParallel(model)
                    model.load_state_dict(checkpoint['model_state_dict'])
                else:
                    state_dict = {k.replace('module.', ''): v for k, v in checkpoint['model_state_dict'].items()}
                    model.load_state_dict(state_dict)     
            else:
                model.load_state_dict(checkpoint['model_state_dict'])
            
            self.logger.info("Model loaded successfully")
            return model
            
        except Exception as e:
            self.logger.error(f"Failed to load model: {str(e)}")
            raise

    def _normalize_features(self, features: torch.Tensor) -> torch.Tensor:
        """
        normalize the features
        
        Args:
            features: input features
            
        Returns:
            torch.Tensor: normalized features
        """
        return F.normalize(features.float(), p=2, dim=1)

    def load_predict_data(self, 
                         img_features_path: str,
                         freq_features_path: str,
                         batch_size: int = 256,
                         num_workers: int = 0) -> DataLoader:
        """
        加载预测数据
        
        Args:
            img_features_path: image feature file path
            freq_features_path: frequency feature file path
            batch_size: batch size
            num_workers: number of workers for data loader
            
        Returns:
            DataLoader: prediction data loader
        """
        try:
            self.logger.info("Loading prediction data...")
            
            # load features
            img_features = torch.load(img_features_path)
            img_features = self._normalize_features(img_features)
            
            freq_features = torch.load(freq_features_path)
            freq_features = self._normalize_features(freq_features)
            
            # create dataset and loader
            predict_dataset = TensorDataset(img_features, freq_features)
            predict_loader = DataLoader(
                predict_dataset,
                batch_size=batch_size,
                shuffle=False,
                num_workers=num_workers
            )
            
            self.logger.info(f"Loaded {len(predict_dataset)} samples for prediction")
            return predict_loader
            
        except Exception as e:
            self.logger.error(f"Failed to load prediction data: {str(e)}")
            raise

    def predict(self, predict_loader: DataLoader) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]:
        """
        predict the labels
        
        Args:
            predict_loader: prediction data loader
            return_probs: whether to return the probabilities
            
        Returns:
            Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]: predicted labels or (labels, probabilities)
        """
        try:
            self.logger.info("Starting prediction...")
            predictions = []
            probabilities = []
            
            with torch.no_grad():
                for batch_x, batch_x_kmer in tqdm(predict_loader, desc="Predicting"):
                    # move data to device
                    batch_x = batch_x.to(self.device)
                    batch_x_kmer = batch_x_kmer.to(self.device)
                    
                    # predict
                    outputs = self.model(batch_x, batch_x_kmer)
                    probs = F.softmax(outputs, dim=1)
                    predicted = probs.argmax(dim=1)
                    
                    # save results
                    predictions.extend(predicted.cpu().numpy())
                    probabilities.extend(probs.cpu().numpy())
            
            predictions = np.array(predictions)
            self.logger.info("Prediction completed")
            
            probabilities = np.array(probabilities)
            return predictions, probabilities
            
        except Exception as e:
            self.logger.error(f"Prediction failed: {str(e)}")
            raise

    def predict_and_save(self,
                        predict_loader: DataLoader,
                        output_dir: str,
                        seq_names: str):
        """
        predict and save the results
        
        Args:
            predict_loader: prediction data loader
            output_dir: output directory
        """
        try:
            if not os.path.exists(output_dir):
                os.makedirs(output_dir, exist_ok=True)
                self.logger.info(f"Output folder -- {output_dir} created successfully")

            predictions, probabilities = self.predict(predict_loader)
            predictions = np.where(probabilities[:, 1] > self.threshold, 1, 0)
            seq_names = pd.read_csv(seq_names, header=None).values.flatten()

            df = pd.DataFrame({'seq_name': seq_names, 'preds': predictions})

            df.to_csv(os.path.join(output_dir, 'is_LTR_deep.txt'), index=False, sep='\t', header=False)

            pd.DataFrame({'seq_name': seq_names, 'preds': predictions}).to_csv(os.path.join(output_dir, 'deep_preds.csv'), index=False, header=False, sep='\t')
            self.logger.info(f"Predictions saved to {output_dir}")

            prob_path = os.path.join(output_dir, 'deep_probs.csv')
            pd.DataFrame({'seq_name': seq_names, 'probs': probabilities[:,1]}).to_csv(prob_path, index=False, header=False, sep='\t')
            self.logger.info(f"Probabilities saved to {prob_path}")
                
        except Exception as e:
            self.logger.error(f"Failed to save predictions: {str(e)}")
            raise

def main():
    parser = argparse.ArgumentParser(description="HybridLTR DeepLearning mudule")
    parser.add_argument("--model_path", type=str, required=True, help="Model file path")
    parser.add_argument("--img_features", type=str, required=True, help="Image feature file path")
    parser.add_argument("--freq_features", type=str, required=True, help="Frequency feature file path")
    parser.add_argument("--seq_names", type=str, required=True, help="Sequence name file path")
    parser.add_argument("--output_dir", type=str, required=True, help="Output file path")
    parser.add_argument("--batch_size", type=int, default=256, help="Batch size")
    parser.add_argument("--threshold", type=float, default=0.9, help="greater than threshold is LTR")
    parser.add_argument("--device", type=str, default='cpu', help="Device type (cuda/cpu)")
    parser.add_argument("--threads", type=int, required=True, help="Number of threads")


    
    args = parser.parse_args()
    
    try:
        # create predictor instance
        predictor = CNNPredictor(args.model_path, args.threshold, args.device)
        
        # load prediction data
        predict_loader = predictor.load_predict_data(
            args.img_features,
            args.freq_features,
            args.batch_size
        )
        
        # predict and save results
        predictor.predict_and_save(
            predict_loader,
            args.output_dir,
            args.seq_names,
        )
        
    except Exception as e:
        logging.error(f"Prediction failed: {str(e)}")
        raise

if __name__ == '__main__':
    main()
