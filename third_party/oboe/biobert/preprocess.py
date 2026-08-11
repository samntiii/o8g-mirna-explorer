#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
from transformers import BertTokenizer
from torch.utils.data import Dataset, DataLoader
import torch

file_path = 'o8G_1.csv'  

def preprocess_data(file_path, tokenizer, max_length=512):
    data = pd.read_csv(file_path)
    
    class RNADataset(Dataset):
        def __init__(self, dataframe, tokenizer, max_length=512):
            self.dataframe = dataframe
            self.tokenizer = tokenizer
            self.max_length = max_length
            self.sequences = dataframe['sequence'].values
            self.labels = dataframe['label'].values

        def __len__(self):
            return len(self.sequences)

        def __getitem__(self, index):
            sequence = self.sequences[index]
            label = self.labels[index]

            encoding = self.tokenizer(
                sequence,
                truncation=True,
                padding='max_length',
                max_length=self.max_length,
                return_tensors='pt'
            )

            return {
                'input_ids': encoding['input_ids'].squeeze(),
                'attention_mask': encoding['attention_mask'].squeeze(),
                'labels': torch.tensor(label, dtype=torch.long)
            }

    return RNADataset(data, tokenizer)

tokenizer = BertTokenizer.from_pretrained('dmis-lab/biobert-base-cased-v1.1')

file_path = 'o8G_1.csv'  
dataset = preprocess_data(file_path, tokenizer)
dataloader = DataLoader(dataset, batch_size=16, shuffle=True)

