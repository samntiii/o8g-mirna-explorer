#!/usr/bin/env python
# coding: utf-8

# In[1]:


import torch
from torch.optim import AdamW
from transformers import BertTokenizer, BertForSequenceClassification, get_linear_schedule_with_warmup
from torch.utils.data import Dataset, DataLoader
from sklearn.metrics import accuracy_score
import pandas as pd
from torch.optim.lr_scheduler import StepLR

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

def train_model(dataloader):
    model = BertForSequenceClassification.from_pretrained('bert-base-uncased', num_labels=2)
    device = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
    model.to(device)

    optimizer = AdamW(model.parameters(), lr=2e-5)
    scheduler = StepLR(optimizer, step_size=3, gamma=0.1)

    epochs = 3
    for epoch in range(epochs):
        model.train()
        running_loss = 0.0
        
        for batch in dataloader:
            input_ids = batch['input_ids'].to(device)
            attention_mask = batch['attention_mask'].to(device)
            labels = batch['labels'].to(device)

            optimizer.zero_grad()

            outputs = model(input_ids, attention_mask=attention_mask, labels=labels)
            loss = outputs.loss

            loss.backward()
            optimizer.step()

            running_loss += loss.item()

        print(f"Epoch {epoch+1}/{epochs} - Loss: {running_loss/len(dataloader)}")
        scheduler.step()

    model.save_pretrained('./finetuned_bert_1')

if __name__ == "__main__":

    tokenizer = BertTokenizer.from_pretrained('bert-base-uncased')

    file_path = 'o8G_1.csv'  

    dataset = preprocess_data(file_path, tokenizer)
    dataloader = DataLoader(dataset, batch_size=16, shuffle=True)

    train_model(dataloader)


# In[ ]:




