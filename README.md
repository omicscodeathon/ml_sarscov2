# Machine learning approach for classification of SARS-CoV-2 variants

### Team Members
1. [Mike Mwanga](https://github.com/mikemwanga)
2. [Evans Mudibo](https://github.com/mudiboevans)
3. [Hesbon Omwandho](https://github.com/hesbornomwandho)
4. [Olaitan I. Awe](https://github.com/laitanawe)
5. [Bonface Onyango](https://github.com/bonfaceonyango)

## Background <br>
The emergence and rapid spread of coronavirus disease 2019 (COVID-19) caused by severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2)is a major and urgent threat to global health. Upon infection, patients exhibit several significant symptoms including fever, cough, diarrhea which are more severe in adults with chronic illness, and shortness of breath[1].

Currently, PCR and sequencing methods have been deployed to identify and classify SARS-CoV-2 variants[2]. However, the design of primers and probes limits these methods to identify only known variants and require continuous update to cater for new species or strains that are clinically relevant.  High-throughput sequencing provides an unbiased identification of virus molecules present in samples and can be used for systematic classification of SARS-CoV-2 variants. However, this approach requires large scaled reference databases of variants to compare against resulting to considerable computing requirements.

Non-clinical techniques such as machine learning (ML) have shown potential for diagnosing, containing and therapeutic monitoring of many diseases. ML based tools can identify and extract the important sequence features for sequence classification in a computationally efficient manner. Such, have been successfully used previously in solving sequence classification problems and some have been applied in identification and classification of SARS-CoV-2. Accurate classification will eventually reduce the huge burden on the limited health system while providing efficient and timely monitoring of transmission and spread patters of SARS-CoV-2.In this work we intend to employ Convolutional Neural Network machine learning architecture on publicly genomics dataset to classify SARS-CoV-2 variants. CNN models are an advancement of machine learning which are technically advanced due to their automated feature extraction.These methods produce excellent results using the concepts of Natural Language Processing(NLP). CNN has been used extensively for DNA classification,protein structure prediction and specifically viral classification.

## Aims and Objectives <br>
To develop a CNN machine learning model that accurately classifies SARS-CoV-2 variants

## Significance of Study <br>
COVID-19 has brought with its immense burden to the healthcare system globally. Currently, vaccination is the only available control measure for control and spread of COVID-19. Accurate classification will enable effective monitoring and tracking of SARS-CoV-2 variants and improve current control and management of the pandemic.

## Methods

### Workflow

Below is a workflow of the methods used in this study <br>

![image](https://github.com/omicscodeathon/ml_sarscov2/blob/main/figures/workflow.png)

### Dataset
Classification of SARS-Cov-2 variants will make use of publicly available representative spike-protein sequence data downloaded from the GISADI Database.

### Data preprocessing of dataset
Python libraries will be used in data processing and exploratory anaysis before feeding to the model. The input features will be the DNA sequence, while the output feature will be the class of SARS-CoV-2 variant. One hot encoding will be used to transform the output feature values into a binary matrix. The dataset will be split into training and testing datasets.

![one_hot](https://github.com/omicscodeathon/ml_sarscov2/blob/main/figures/one_hot_encoder.png)

ML algorithms do make use of numerical data and since DNA sequences are in categorical form, some form of conversion using readly available tools will be implemented. In this study, one-hot encoding is used to convert the k-merized sequence data into numerical form.

### Classification model
The resulting dataset is fed into a CNN model for feature   extraction and classification. CNN uses convolutionary layers to automatically extract features from a dataset as opposed to other models which require the user to manually extract important features. CNN contains several layers, one input, several hidden and one output layer and each layer contains several neurons and each neuron contains several parameters.




## References <br>
[1]	H. A. Rothan and S. N. Byrareddy, “The epidemiology and pathogenesis of coronavirus disease (COVID-19) outbreak,” Journal of Autoimmunity, vol. 109. 2020. doi: 10.1016/j.jaut.2020.102433.

[2]	M. A. Shereen, S. Khan, A. Kazmi, N. Bashir, and R. Siddique, “COVID-19 infection: Origin, transmission, and characteristics of human coronaviruses,” Journal of Advanced Research, vol. 24. 2020. doi: 10.1016/j.jare.2020.03.005.

[3]	S. Sethi and T. Chakraborty, “Molecular (real-time reverse transcription polymerase chain reaction) diagnosis of SARS-CoV-2 infections: Complexity and challenges,” Journal of Laboratory Medicine, vol. 45, no. 3. 2021. doi: 10.1515/labmed-2020-0135.

[4]	D. A. Mendels et al., “Using artificial intelligence to improve COVID-19 rapid diagnostic test result interpretation,” Proceedings of the National Academy of Sciences of the United States of America, vol. 118, no. 12, Mar. 2021, doi: 10.1073/pnas.2019893118.

[5]	O. P. Singh, M. Vallejo, I. M. El-Badawy, A. Aysha, J. Madhanagopal, and A. A. Mohd Faudzi, “Classification of SARS-CoV-2 and non-SARS-CoV-2 using machine learning algorithms,” Computers in Biology and Medicine, vol. 136, 2021, doi: 10.1016/j.compbiomed.2021.104650.

[6]	G. S. Randhawa, M. P. M. Soltysiak, H. el Roz, C. P. E. de Souza, K. A. Hill, and L. Kari, “Machine learning using intrinsic genomic signatures for rapid classification of novel pathogens: COVID-19 case study,” PLoS ONE, vol. 15, no. 4, Apr. 2020, doi: 10.1371/journal.pone.0232391.

[7]	Z. Bzhalava, A. Tampuu, P. Bała, R. Vicente, and J. Dillner, “Machine Learning for detection of viral sequences in human metagenomic datasets,” BMC Bioinformatics, vol. 19, no. 1, Sep. 2018, doi: 10.1186/s12859-018-2340-x.

[8] C.M. Dasari,R. Bhukya ,"Explainable deep neural networks for novel viral genome prediction", vol. 52.2021,doi: 10.1007/s10489-021-02572-3
