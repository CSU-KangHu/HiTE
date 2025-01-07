import argparse
from hybridLTR_feature_extractor import FeatureExtractor
from hybridLTR_predicter import CNNPredictor
import logging
import os
import pandas as pd
import torch

def extractor(args):
    try:
        # create feature extractor instance
        feature_extractor = FeatureExtractor(args.matrix_dir, args.threads)

        if not os.path.exists(args.feature_output_dir):
            os.makedirs(args.feature_output_dir, exist_ok=True)
            feature_extractor.logger.info(f"Feature output folder -- {args.feature_output_dir} created successfully")

        # extract features and save results
        img_features, freq_features, seq_names = feature_extractor.extract_main()
        torch.save(img_features, os.path.join(args.feature_output_dir, 'img_features.pt'))
        torch.save(freq_features, os.path.join(args.feature_output_dir, 'freq_features.pt'))
        pd.DataFrame(seq_names).to_csv(os.path.join(args.feature_output_dir, 'seq_names.txt'), index=False, header=False)

        feature_extractor.cleanup()
    except Exception as e:
        print(f"Feature extraction failed: {str(e)}")
        raise

def predictor(args):
    try:
        # create predictor instance
        predictor = CNNPredictor(args.model_path, args.threshold, args.device)
        
        # Load prediction data
        predict_loader = predictor.load_predict_data(
            args.img_features,
            args.freq_features,
            args.batch_size
        )
        
        # Predict and save results
        predictor.predict_and_save(
            predict_loader,
            args.output_dir,
            args.seq_names,
        )
        
    except Exception as e:
        logging.error(f"Prediction failed: {str(e)}")
        raise

def main():
    parser = argparse.ArgumentParser(description="HybridLTR DeepLearning mudule")
    parser.add_argument("--matrix_dir", type=str, required=True, help="The path of the flank alignment matrix")
    parser.add_argument("--threads", type=int, required=True, help="Number of threads")
    parser.add_argument("--feature_output_dir", type=str, required=True, help="Feature output file path")

    parser.add_argument("--model_path", type=str, required=True, help="Model file path")
    parser.add_argument("--img_features", type=str, required=True, help="Image feature file path")
    parser.add_argument("--freq_features", type=str, required=True, help="Frequency feature file path")
    parser.add_argument("--seq_names", type=str, required=True, help="Sequence name file path")
    parser.add_argument("--output_dir", type=str, required=True, help="Output file path")
    parser.add_argument("--batch_size", type=int, default=256, help="Batch size")
    parser.add_argument("--threshold", type=float, default=0.9, help="greater than threshold is LTR")
    parser.add_argument("--device", type=str, default='cpu', help="Device type (cuda/cpu)")
        
    args = parser.parse_args()

    extractor(args)

    predictor(args)

if __name__ == '__main__':
    main()