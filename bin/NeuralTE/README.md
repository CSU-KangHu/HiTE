# NeuralTE
[![GitHub](https://img.shields.io/badge/python-3-blue)](https://www.python.org/)
[![GitHub](https://img.shields.io/badge/license-GPL--3.0-green)](https://github.com/CSU-KangHu/NeuralTE/blob/master/LICENSE)
[![Conda](https://img.shields.io/badge/Conda-support-yellow)](https://docs.conda.io/en/latest/)

`NeuralTE` uses a Convolutional Neural Network (CNN) to classify transposable elements (TEs) at the **superfamily** level, based on the **sequence** and **structural features** of transposons.

It is recommended that the TE library to be classified consists of full-length TEs. Many TE libraries often  divide LTR retrotransposons into terminal and internal sequences, such as `Copia-62_PHord-LTR` and `Copia-62_PHord-I` in Repbase, and we suggest restoring them to full-length LTRs before classification. 

We recommend using [HiTE](https://github.com/CSU-KangHu/HiTE) for generating full-length TE libraries, as it has the capability to identify more full-length TEs, including novel ones.
For fragmented TEs, you can use [RepeatClassifier](https://github.com/Dfam-consortium/RepeatModeler/blob/master/RepeatClassifier) configured with a [complete Dfam library](https://www.repeatmasker.org/RepeatMasker/) for classification using homology searches.

## Table of Contents
- [Installation](#install)
- [Pre-trained models](#models)
- [Specify the GPUs](#gpu)
- [Demo data](#demo)
- [Train a new model](#train_model)
- [Experiment reproduction](#reproduction)
- [Usage](#cmd)

## <a name="install"></a>Installation
NeuralTE is built on [Python3](https://www.python.org/) and [Keras](https://keras.io/).
   - Prerequisites: \
       [Python3](https://www.python.org/) (version=3.8.16)\
       [CUDA Toolkit](https://anaconda.org/anaconda/cudatoolkit) (version>=11.2, for GPU only)
   - Dependencies: \
       [tensorflow](https://www.tensorflow.org/) (version=2.6.0) \
       [Keras](https://keras.io/) (version=2.6.0) \
       [numpy](http://www.numpy.org/) \
       [biopython](https://biopython.org/) \
       [scikit-learn](https://scikit-learn.org/stable/) \
       [matplotlib](https://matplotlib.org/) \
       [seaborn](https://seaborn.pydata.org/) \
       [rmblast](https://www.repeatmasker.org/rmblast/) \
       [openpyxl](https://openpyxl.readthedocs.io/)

#### System Requirements
`NeuralTE` requires a standard computer to use the Convolutional Neural Network (CNN). Using GPU could acceralate the process of TE classification.

Recommended Hardware requirements: 40 CPU processors, 128 GB RAM.

Recommended OS: (Ubuntu 16.04, CentOS 7, etc.)

```sh
git clone https://github.com/CSU-KangHu/NeuralTE.git ## Alternatively, you can download the zip file directly from the repository.
cd NeuralTE
chmod +x tools/*

conda install mamba -c conda-forge
mamba env create --name NeuralTE -f environment.yml
conda activate NeuralTE
```

## <a name="models"></a>Pre-trained models
See [models](/models):
* [NeuralTE model](/models/NeuralTE_model.h5): This model employs features such as k-mer occurrences, terminals, TE domain, and 5bp terminals, trained using Repbase version 28.06.

* [NeuralTE-TSDs model](/models/NeuralTE-TSDs_model.h5): This model incorporates features like k-mer occurrences, terminals, TE domain, 5bp terminals, and **target site duplications (TSDs)**. It was trained using partial species data (493 species) from Repbase version 28.06. Please note that this model should be used in conjunction with the corresponding genome assembly of the species.

## <a name="gpu"></a>Specify the GPUs (Skipping when using CPUs)
```sh
  --start_gpu_num start_gpu_num
                        The starting index for using GPUs. default = [ 0 ]
  --use_gpu_num use_gpu_num
                        Specifying the number of GPUs in use. default = [ 1 ]
```
For example, `--start_gpu_num 0` and `--use_gpu_num 2` indicates a total of two GPUs to be used, with the assigned GPU indices being `gpu:0` and `gpu:1`. 

Default configurations are set in [gpu_config.py](/configs/gpu_config.py).

## <a name="demo"></a>Demo data

Please refer to [demo](/demo) for some demo data to play with:
* _test.fa_: demo TE library.
* _genome.fa_: demo genome sequence.

```sh
# 1.Classify TE library without genome
# Inputs: 
#       --data: TE library to be classified.
#       --model_path: Pre-trained NeuralTE model without using TSDs features.
#       --out_dir: Output directory. The `--out_dir` should not be the same as the directory 
#                 where the `--data` file is located.
# Outputs: 
#       classified.info: Classification labels corresponding to TE names.
#       classified_TE.fa: Classified TE library.
python ${pathTo}/NeuralTE/src/Classifier.py \
 --data ${pathTo}/NeuralTE/demo/test.fa \
 --model_path ${pathTo}/NeuralTE/models/NeuralTE_model.h5 \
 --out_dir ${out_dir} \
 --thread ${threads_num}
 # e.g., my command: 
 # python /home/hukang/NeuralTE/src/Classifier.py \
 # --data /home/hukang/NeuralTE/demo/test.fa \
 # --model_path /home/hukang/NeuralTE/models/NeuralTE_model.h5 \
 # --out_dir /home/hukang/NeuralTE/demo/work \
 # --thread 40
 
 
 # 2.Classify the TE library with genome
 #       test.fa: TE library to be classified 
 #       genome.fa: The genome assembly corresponding to TE library
 #       NeuralTE-TSDs_model.h5: Pre-trained Neural TE model using TSDs features
 # outputs: 
 #       classified.info: Classification labels corresponding to TE names
 #       classified_TE.fa: Classified TE library
python ${pathTo}/NeuralTE/src/Classifier.py \
 --data ${pathTo}/NeuralTE/demo/test.fa \
 --genome ${pathTo}/NeuralTE/demo/genome.fa \
 --use_TSD 1 \
 --model_path ${pathTo}/NeuralTE/models/NeuralTE-TSDs_model.h5 \
 --out_dir ${out_dir} \
 --thread ${threads_num}
 # e.g., my command: 
 # python /home/hukang/NeuralTE/src/Classifier.py \
 # --data /home/hukang/NeuralTE/demo/test.fa \
 # --genome /home/hukang/NeuralTE/demo/genome.fa \
 # --use_TSD 1 \
 # --model_path /home/hukang/NeuralTE/models/NeuralTE-TSDs_model.h5 \
 # --out_dir /home/hukang/NeuralTE/demo/work \
 # --thread 40
```

## <a name="train_model"></a>Train a new model
- Prerequisites: \
       [Repbase*.fasta.tar.gz](https://www.girinst.org/server/RepBase/index.php) (version>=28.06)\
       [Genomes](https://www.ncbi.nlm.nih.gov/) (Optional, for TSDs model only)
- For more detailed information, please refer to [Experiment reproduction](#reproduction).

```sh
# 0. Preprocess Repbase database (including merging subfiles, concatenating LTR terminal and internal sequences, filtering incomplete LTR transposons, etc.)
# Inputs:
#        repbase_dir: Directory containing all Repbase files
#        out_dir: Output directory containing preprocessed results
# Outputs:
#        all_repbase.ref: Merged sequence of all Repbase databases
python ${pathTo}/utils/preprocess_repbase.py \
 --repbase_dir ${pathTo}/RepBase${version}.fasta \
 --out_dir ${out_dir}
 # e.g., my command: 
 # python /home/hukang/NeuralTE/utils/preprocess_repbase.py \
 # --repbase_dir /home/hukang/RepBase28.06.fasta/ \
 # --out_dir /home/hukang/test/
 
 
# 1. Splitting train and test datasets
# Inputs:
#        data_path: All Repbase database sequences
#        out_dir: Output directory after dataset partition
#        ratio: Ratio of training set to test set
# Outputs:
#        train.ref: 80% of all Repbase database sequences for training
#        test.ref: 20% of all Repbase database sequences for testing
python ${pathTo}/utils/split_train_test.py \
 --data_path ${Step0_out_dir}/all_repbase.ref \
 --out_dir ${out_dir} \
 --ratio 0.8
 # e.g., my command: 
 # python /home/hukang/NeuralTE/utils/split_train_test.py \
 # --data_path /home/hukang/test/all_repbase.ref \
 # --out_dir /home/hukang/test/ \
 # --ratio 0.8


 # 2.Train a new NeuralTE Model
 # Inputs: 
 #        train.ref: training set
 # Outputs: 
 #        model_${time}.h5: Generate h5 format file in the ${pathTo}/NeuralTE/models directory
python ${pathTo}/NeuralTE/src/Trainer.py \
 --data ${Step1_out_dir}/train.ref \
 --is_train 1 \
 --is_predict 0 \
 --use_TSD 0 \
 --out_dir ${out_dir} \
 --thread ${threads_num} \
 --start_gpu_num ${start_gpu_num} \
 --use_gpu_num ${use_gpu_num}
 # e.g., my command: 
 # python /home/hukang/NeuralTE/src/Trainer.py \
 # --data /home/hukang/NeuralTE/data/train.ref \
 # --is_train 1 \
 # --is_predict 0 \
 # --use_TSD 0 \
 # --out_dir /home/hukang/NeuralTE/work \
 # --thread 40
 # --start_gpu_num 0 \
 # --use_gpu_num 1
 
 # Replace original model
 cd ${pathTo}/NeuralTE/models && mv model_${time}.h5 NeuralTE_model.h5
 
 
 # 3.Train a new NeuralTE-TSDs Model
 # Inputs: 
 #        train.ref: training set
 #        genome.info: Modify the 'genome.info' file in the ${pathTo}/NeuralTE/data directory. Ensure that 'Scientific Name' corresponds to the species names in `train.ref`, and 'Genome Path' should be an absolute path.
 # Outputs: 
 #        model_${time}.h5: Generate h5 format file in the ${pathTo}/NeuralTE/models directory
python ${pathTo}/NeuralTE/src/Trainer.py \
 --data ${Step1_out_dir}/train.ref \
 --genome ${pathTo}/NeuralTE/data/genome.info \
 --is_train 1 \
 --is_predict 0 \
 --use_TSD 1 \
 --out_dir ${out_dir} \
 --thread ${threads_num} \
 --start_gpu_num ${start_gpu_num} \
 --use_gpu_num ${use_gpu_num}
 # e.g., my command: 
 # python /home/hukang/NeuralTE/src/Trainer.py \
 # --data /home/hukang/NeuralTE/data/train.ref \
 # --genome /home/hukang/NeuralTE/data/genome.info \
 # --is_train 1 \
 # --is_predict 0 \
 # --use_TSD 1 \
 # --out_dir /home/hukang/NeuralTE/work \
 # --thread 40 \
 # --start_gpu_num 0 \
 # --use_gpu_num 1
 
 # Replace original model
cd ${pathTo}/NeuralTE/models && mv model_${time}.h5 NeuralTE-TSDs_model.h5
```

<!--
## <a name="libraries"></a>TE libraries generated by HiTE, EDTA, and RepeatModeler2 

To reduce redundancy, many tools partition LTR retrotransposons into terminal and internal sequences. We recommend restoring them to full-length LTRs before classification.
In the [demo](/demo) directory, authentic TE libraries for HiTE, EDTA, and RepeatModeler2 are provided. You can attempt to classify them using NeuralTE following the instructions below.
* For the HiTE output, you can directly use the file `confident_TE.cons.fa` as input for NeuralTE. NeuralTE will automatically concatenate the terminal and internal sequences of LTR retrotransposons.
```sh
python ${pathTo}/NeuralTE/src/Classifier.py \
 --data ${pathTo}/NeuralTE/demo/confident_TE.cons.fa \
 --model_path ${pathTo}/NeuralTE/models/NeuralTE_model.h5 \
 --out_dir ${out_dir} \
 --thread ${threads_num}
 # e.g., my command: python /home/hukang/NeuralTE/src/Classifier.py \
 # --data /home/hukang/NeuralTE/demo/confident_TE.cons.fa \
 # --model_path /home/hukang/NeuralTE/models/NeuralTE_model.h5 \
 # --out_dir /home/hukang/NeuralTE/demo/work \
 # --thread 40
```

* For the EDTA output, we suggest utilizing the complete TE library, `$genome.mod.EDTA.final/$genome.mod.EDTA.intact.fa`, as input for NeuralTE.
```sh
python ${pathTo}/NeuralTE/src/Classifier.py \
 --data ${pathTo}/NeuralTE/demo/GCF_000001735.4_TAIR10.1_genomic.rename.fna.mod.EDTA.intact.fa \
 --model_path ${pathTo}/NeuralTE/models/NeuralTE_model.h5 \
 --out_dirout_dir ${out_dir} \
 --thread ${threads_num}
 # e.g., my command: python /home/hukang/NeuralTE/src/Classifier.py \
 # --data /home/hukang/NeuralTE/demo/GCF_000001735.4_TAIR10.1_genomic.rename.fna.mod.EDTA.intact.fa \
 # --model_path /home/hukang/NeuralTE/models/NeuralTE_model.h5 \
 # --out_dir /home/hukang/NeuralTE/demo/work \
 # --thread 40
```

* For the RepeatModeler2 output, we advise using the `NeuralTE/utils/reName_RM2.py` script to rename the headers before utilizing it as input for NeuralTE.
```sh
# 1. Rename the LTR retrotransposons in the TE library output from RepeatModeler2, with the output header format: `RM2_intactLTR_114-LTR` and `RM2_intactLTR_114-I`
python ${pathTo}/NeuralTE/utils/reName_RM2.py \
 -i ${pathTo}/${species}-families.fa \
 -o ${pathTo}/${species}-families.rename.fa
  # e.g., my command: python /home/hukang/NeuralTE/utils/reName_RM2.py \
 # -i /home/hukang/NeuralTE/demo/rice-families.fa \
 # -o /home/hukang/NeuralTE/demo/rice-families.rename.fa
 
 # 2. Classify RepeatModeler2 library
 python ${pathTo}/NeuralTE/src/Classifier.py \
 --data ${pathTo}/${species}-families.rename.fa \
 --model_path ${pathTo}/NeuralTE/models/NeuralTE_model.h5 \
 --out_dir ${out_dir} \
 --thread ${threads_num}
 # e.g., my command: python /home/hukang/NeuralTE/src/Classifier.py \
 # --data /home/hukang/NeuralTE/demo/rice-families.rename.fa \
 # --model_path /home/hukang/NeuralTE/models/NeuralTE_model.h5 \
 # --out_dir /home/hukang/NeuralTE/demo/work \
 # --thread 40
```
-->

## <a name="reproduction"></a>Experiment reproduction

All experimental results in the manuscript of NeuralTE can be reproduced step by step through [Experiment reproduction](https://github.com/CSU-KangHu/NeuralTE/wiki/Experiment-reproduction).

## <a name="cmd"></a>Usage
#### 1. Classify TE library
```shell
usage: Classifier.py [-h] --data data --out_dir output_dir [--use_TSD use_TSD] [--is_predict is_predict] [--start_gpu_num start_gpu_num] [--use_gpu_num use_gpu_num] [--keep_raw keep_raw] [--genome genome] [--species species] [--model_path model_path]
                     [--use_kmers use_kmers] [--use_terminal use_terminal] [--use_minority use_minority] [--use_domain use_domain] [--use_ends use_ends] [--is_wicker is_wicker] [--is_plant is_plant] [--threads thread_num] [--internal_kmer_sizes internal_kmer_sizes]
                     [--terminal_kmer_sizes terminal_kmer_sizes]

########################## NeuralTE, version 1.0.0 ##########################

optional arguments:
  -h, --help            show this help message and exit
  --data data           Input fasta file used to predict, header format: seq_name label species_name, refer to "data/test.example.fa" for example.
  --out_dir output_dir  Output directory, store temporary files
  --use_TSD use_TSD     Whether to use TSD features, 1: true, 0: false. default = [ 0 ]
  --is_predict is_predict
                        Enable prediction mode, 1: true, 0: false. default = [ 1 ]
  --start_gpu_num start_gpu_num
                        The starting index for using GPUs. default = [ 0 ]
  --use_gpu_num use_gpu_num
                        Specifying the number of GPUs in use. default = [ 1 ]
  --keep_raw keep_raw   Whether to retain the raw input sequence, 1: true, 0: false; only save species having TSDs. default = [ 0 ]
  --genome genome       Genome path, use to search for TSDs
  --species species     Which species does the TE library to be classified come from.
  --model_path model_path
                        Input the path of trained model, absolute path.
  --use_kmers use_kmers
                        Whether to use kmers features, 1: true, 0: false. default = [ 1 ]
  --use_terminal use_terminal
                        Whether to use LTR, TIR terminal features, 1: true, 0: false. default = [ 1 ]
  --use_minority use_minority
                        Whether to use minority features, 1: true, 0: false. default = [ 0 ]
  --use_domain use_domain
                        Whether to use domain features, 1: true, 0: false. default = [ 1 ]
  --use_ends use_ends   Whether to use 5-bp terminal ends features, 1: true, 0: false. default = [ 1 ]
  --is_wicker is_wicker
                        Use Wicker or RepeatMasker classification labels, 1: Wicker, 0: RepeatMasker. default = [ 1 ]
  --is_plant is_plant   Is the input genome of a plant? 0 represents non-plant, while 1 represents plant. default = [ 0 ]
  --threads thread_num  Input thread num, default = [ 104 ]
  --internal_kmer_sizes internal_kmer_sizes
                        The k-mer size used to convert internal sequences to k-mer frequency features, default = [ [5] MB ]
  --terminal_kmer_sizes terminal_kmer_sizes
                        The k-mer size used to convert terminal sequences to k-mer frequency features, default = [ [3, 4] ]
```
#### 2. Train a new model
```shell
usage: Trainer.py [-h] --data data --out_dir output_dir --use_TSD use_TSD --is_train is_train --is_predict is_predict [--start_gpu_num start_gpu_num] [--use_gpu_num use_gpu_num] [--only_preprocess only_preprocess] [--keep_raw keep_raw] [--genome genome]
                  [--use_kmers use_kmers] [--use_terminal use_terminal] [--use_minority use_minority] [--use_domain use_domain] [--use_ends use_ends] [--threads thread_num] [--internal_kmer_sizes internal_kmer_sizes] [--terminal_kmer_sizes terminal_kmer_sizes]
                  [--cnn_num_convs cnn_num_convs] [--cnn_filters_array cnn_filters_array] [--cnn_kernel_sizes_array cnn_kernel_sizes_array] [--cnn_dropout cnn_dropout] [--batch_size batch_size] [--epochs epochs] [--use_checkpoint use_checkpoint]

########################## NeuralTE, version 1.0.0 ##########################

optional arguments:
  -h, --help            show this help message and exit
  --data data           Input fasta file used to train model, header format: seq_name label species_name, refer to "data/train.example.fa" for example.
  --out_dir output_dir   Output directory, store temporary files
  --use_TSD use_TSD     Whether to use TSD features, 1: true, 0: false. default = [ 0 ]
  --is_train is_train   Enable train mode, 1: true, 0: false. default = [ 0 ]
  --is_predict is_predict
                        Enable prediction mode, 1: true, 0: false. default = [ 1 ]
  --start_gpu_num start_gpu_num
                        The starting index for using GPUs. default = [ 0 ]
  --use_gpu_num use_gpu_num
                        Specifying the number of GPUs in use. default = [ 1 ]
  --only_preprocess only_preprocess
                        Whether to only perform data preprocessing, 1: true, 0: false.
  --keep_raw keep_raw   Whether to retain the raw input sequence, 1: true, 0: false; only save species having TSDs. default = [ 0 ]
  --genome genome       Genome path, use to search for TSDs
  --use_kmers use_kmers
                        Whether to use kmers features, 1: true, 0: false. default = [ 1 ]
  --use_terminal use_terminal
                        Whether to use LTR, TIR terminal features, 1: true, 0: false. default = [ 1 ]
  --use_minority use_minority
                        Whether to use minority features, 1: true, 0: false. default = [ 0 ]
  --use_domain use_domain
                        Whether to use domain features, 1: true, 0: false. default = [ 1 ]
  --use_ends use_ends   Whether to use 5-bp terminal ends features, 1: true, 0: false. default = [ 1 ]
  --threads thread_num  Input thread num, default = [ 104 ]
  --internal_kmer_sizes internal_kmer_sizes
                        The k-mer size used to convert internal sequences to k-mer frequency features, default = [ [5] MB ]
  --terminal_kmer_sizes terminal_kmer_sizes
                        The k-mer size used to convert terminal sequences to k-mer frequency features, default = [ [3, 4] ]
  --cnn_num_convs cnn_num_convs
                        The number of CNN convolutional layers. default = [ 3 ]
  --cnn_filters_array cnn_filters_array
                        The number of filters in each CNN convolutional layer. default = [ [16, 16, 16] ]
  --cnn_kernel_sizes_array cnn_kernel_sizes_array
                        The kernel size in each of CNN convolutional layer. default = [ [7, 7, 7] ]
  --cnn_dropout cnn_dropout
                        The threshold of CNN Dropout. default = [ 0.5 ]
  --batch_size batch_size
                        The batch size in training model. default = [ 32 ]
  --epochs epochs       The number of epochs in training model. default = [ 50 ]
  --use_checkpoint use_checkpoint
                        Whether to use breakpoint training. 1: true, 0: false. The model will continue training from the last failed parameters to avoid training from head. default = [ 0 ]
```


