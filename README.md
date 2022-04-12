# Machine learning approach for classification of SARS-CoV-2 variants

### Team Members
1. [Mike Mwanga](https://github.com/mikemwanga)
2. [Evans Mudibo](https://github.com/mudiboevans)
3. [Hesbon Omwandho](https://github.com/hesbornomwandho)
4. [Olaitan I. Awe](https://github.com/laitanawe)
5. [Bonface Onyango](https://github.com/bonfaceonyango)

## Background <br>
The emergence and rapid spread of coronavirus disease 2019 (COVID-19) caused by severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2) as a potentially fatal disease is a major and urgent threat to global health. COVID-19 is transmitted through direct contact with an infected person via sneezing and coughing and has no medically approved vaccine or medication[1], [2]. Upon infection, patients exhibit several significant symptoms including fever, cough, diarrhea which are more severe in adults with chronic illness, and shortness of breath[2].

Currently, PCR and sequencing methods have been deployed to identify and classify SARS-CoV-2 variants[3]. However, the design of primers and probes limits these methods to identify only known variants and require continuous update to cater for new species or strains that are clinically relevant.  High-throughput sequencing provides an unbiased identification of virus molecules present in samples and can be used for systematic classification of SARS-CoV-2 variants. However, this approach requires large scaled reference databases of variants to compare against resulting to considerable computing requirements.

Non-clinical techniques such as machine learning (ML), and other artificial intelligence (AI) techniques can be implemented in diagnosis, classification and eventually control of COVID-19[4], [5]. ML based tools can identify and extract the important sequence features for sequence classification in a computationally efficient manner. Such, have been successfully used previously in solving sequence classification problems[5]. Despite these advances, none has been applied in classifying SARS-CoV-2 variants. Accurate classification of SARS-CoV-2 variants will eventually reduce the huge burden on the limited health system while providing efficient and diagnostic and predictable methods for SARS-CoV-2.
Machine learning (ML) is one of the most advanced AI techniques which has shown potential for diagnosing, containing and therapeutic monitoring of many diseases. In this work we intend to employ machine learning techniques on publicly available epidemiological and genomics dataset to predict infection by SARS-CoV-2 and classify variants.
## Current techniques for viral classification
Various current techniques for classification  of viral genomes can be categorised into alignment-based or machine learning based. 
Alignment-based Techniques
In an alignment-based approach, viral sequence is determined using BLAST,where the sequences are compared to known available sequences in the databases. The viral sequences are then classified based on similarity index. The most common alignment-based techniques are USEARCH,SCUAL and REGA,HHMERS.
## Limitations of Alignment-based methods
The main disadvantage of using alignment-based methodology is that it involves several initial sequence alignments and hyperparameters. Consequently, they are computationally expensive, especially with large datasets. These approaches do not perform well with divergent regions of genomes. To counteract thes challenges, alignment-free models have been adopted
## Alignment-free approaches
Machine learning methods
Various machine learning models have been used to classify viral genomes.The most common alignment-free alignments are VirSorter,ViralFinder. Virtsorter is a probabilistic model tool used to predict viral sequence with or without a reference sequence.VirFinder is a machine learning model that uses k-mer frequency to classify viral genomes.Another recommended machine learning model tool is ViraPipe which used random forest and artificial neural network for classification.To improve on classification,ViraPipe uses relative synonymous codon usage frequency. The model identified two codons (TCG and GCG) which have been shown to have very strong discriminative ability.These methods are still faced with drawbacks.
## Limitations of machine-learning techniques
One major challenge with machine learning models is their inability to extract important hidden features from basic DNA. Most machine learning models rely solely on the extracted features. These features can be identified manually or automated. While manual feature extraction is laborious and requires one to have a deep understanding of how a model works, automation requires much computer memory which still makes machine learning models computationally expensive. Machine learning algorithm is interpretable if the user can deeply understand and observe the model parameters and how the model makes choices on its own. Explainable models are too complex to comprehend and hence require extra techniques to interpret the model performance,evaluation and genome classification.
To counteract machine learning challenges, deep learning models have been currently adopted.
## Deep learning models.
Deep learning models are an advancement of machine learning which are technically advanced due to their automated feature extraction.These methods produce excellent results using the concepts of Natural Language Processing(NLP).
Deep learning has been used extensively for DNA classification,protein structure prediction and genome sequence analysis. Convolutional Neural Networks have been widely proposed for viral classification.AUC-ROC is the metric used in evaluating the performance of the models. Based on AUC-ROC, CNN and LSTM based models have been proven for excellent performance. ViralMiner  is one of the CNN based models that  has been used to classify viral genomes.The model has achieved pattern of 0.905 and frequency branch of 0.917 AUC-ROC. Another CNN model tool, DeepVirFinder, has been used to classify novel viral genomes  and its performance increases with increase in viral sequence length. It achieved 0.978 AUC-ROC in classification of viral sequences of length 3000 bp and . RNN-ViralSeeker is an LSTM model for classification of viral sequences with AUC-ROC of 0.9175[8]



## Aims and Objectives <br>

To classify SARS-CoV-2 variants using genomic sequence data and machine learning

## Significance of Study <br>
COVID-19 has brought with its immense burden to the healthcare system globally. Currently, vaccination is the only available control measure for control and spread of COVID-19. Accurate classification will enable effective monitoring and tracking of SARS-CoV-2 variants and improve in control and management of the pandemic.

## Methods

### Workflow

Below is a workflow of the methods used in this study <br>

![image](https://github.com/omicscodeathon/ml_sarscov2/blob/main/figures/workflow.png)

### Dataset
Classification of SARS-Cov-2 variants will make use of publicly available representative spike-protein sequence data downloaded from the GISADI and/or GenBank Databases.

### Loading and preprocessing of dataset
The dataset will be loaded using pandas. The input features will be the DNA sequence, while the output feature will be the class of SARS-CoV-2 variant. One hot encoding will be used to transform the output feature values into a binary matrix. The dataset will be split into training and testing datasets.

![one_hot](https://github.com/omicscodeathon/ml_sarscov2/blob/main/figures/one_hot_encoder.png)

### Data processing
ML algorithms do make use of numerical data and since DNA sequences are in categorical form, some form of conversion using readly available tools will be implemented. In this study, label encoding and k-mer encoding techniques are used to convert the sequence data into numerical form. Seqeuences will first be converted into k-mers (which size will be appropriate, 3? In relation to codons?). This will result to k-mer patterns specific for each variant. The label-encoding process will implement LabelEncoder(), where, each k-mer is assigned a numerical value in a sequential manner. The resulting 2D sequence representation numerical matrix will be binarized using  LabelBinarizer()

### Classification model
The resulting dataset is fed into a convolutional neural network model (CNN) for feature   extraction. CNN uses convolutionary layers to automatically extract features from a dataset as opposed to other models which require the user to manually extract important features. CNN contains several layers, one input, several hidden and one output layer and each layer contains several neurons and each neuron contains several parameters.



## References <br>
[1]	H. A. Rothan and S. N. Byrareddy, “The epidemiology and pathogenesis of coronavirus disease (COVID-19) outbreak,” Journal of Autoimmunity, vol. 109. 2020. doi: 10.1016/j.jaut.2020.102433.

[2]	M. A. Shereen, S. Khan, A. Kazmi, N. Bashir, and R. Siddique, “COVID-19 infection: Origin, transmission, and characteristics of human coronaviruses,” Journal of Advanced Research, vol. 24. 2020. doi: 10.1016/j.jare.2020.03.005.

[3]	S. Sethi and T. Chakraborty, “Molecular (real-time reverse transcription polymerase chain reaction) diagnosis of SARS-CoV-2 infections: Complexity and challenges,” Journal of Laboratory Medicine, vol. 45, no. 3. 2021. doi: 10.1515/labmed-2020-0135.

[4]	D. A. Mendels et al., “Using artificial intelligence to improve COVID-19 rapid diagnostic test result interpretation,” Proceedings of the National Academy of Sciences of the United States of America, vol. 118, no. 12, Mar. 2021, doi: 10.1073/pnas.2019893118.

[5]	O. P. Singh, M. Vallejo, I. M. El-Badawy, A. Aysha, J. Madhanagopal, and A. A. Mohd Faudzi, “Classification of SARS-CoV-2 and non-SARS-CoV-2 using machine learning algorithms,” Computers in Biology and Medicine, vol. 136, 2021, doi: 10.1016/j.compbiomed.2021.104650.

[6]	G. S. Randhawa, M. P. M. Soltysiak, H. el Roz, C. P. E. de Souza, K. A. Hill, and L. Kari, “Machine learning using intrinsic genomic signatures for rapid classification of novel pathogens: COVID-19 case study,” PLoS ONE, vol. 15, no. 4, Apr. 2020, doi: 10.1371/journal.pone.0232391.

[7]	Z. Bzhalava, A. Tampuu, P. Bała, R. Vicente, and J. Dillner, “Machine Learning for detection of viral sequences in human metagenomic datasets,” BMC Bioinformatics, vol. 19, no. 1, Sep. 2018, doi: 10.1186/s12859-018-2340-x.

[8] C.M. Dasari,R. Bhukya ,"Explainable deep neural networks for novel viral genome prediction", vol. 52.2021,doi: 10.1007/s10489-021-02572-3
