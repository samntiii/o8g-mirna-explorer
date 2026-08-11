#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import torch
from transformers import BertTokenizer, BertForSequenceClassification
from torch.utils.data import Dataset, DataLoader
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, confusion_matrix, roc_auc_score
import pandas as pd

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

def evaluate_model(model, dataloader, device):
    model.eval()
    predictions = []
    true_labels = []
    probs = [] 

    with torch.no_grad():
        for batch in dataloader:
            input_ids = batch['input_ids'].to(device)
            attention_mask = batch['attention_mask'].to(device)
            labels = batch['labels'].to(device)

            outputs = model(input_ids, attention_mask=attention_mask)
            logits = outputs.logits
            softmax_probs = torch.softmax(logits, dim=-1)[:, 1]  # 🔧 取正类概率
            preds = torch.argmax(logits, dim=-1)

            predictions.extend(preds.cpu().numpy())
            true_labels.extend(labels.cpu().numpy())
            probs.extend(softmax_probs.cpu().numpy())  # 🔧

    accuracy = accuracy_score(true_labels, predictions)
    precision = precision_score(true_labels, predictions, average='binary')
    recall = recall_score(true_labels, predictions, average='binary')
    f1 = f1_score(true_labels, predictions, average='binary')

    print(f'Accuracy on test set: {accuracy:.4f}')
    print(f'Precision: {precision:.4f}')
    print(f'Recall: {recall:.4f}')
    print(f'F1 Score: {f1:.4f}')

    from sklearn.metrics import roc_curve, auc
    import numpy as np

    fpr, tpr, thresholds = roc_curve(true_labels, probs)
    roc_auc = auc(fpr, tpr)
    print(f"AUC: {roc_auc:.4f}")

    np.save("fpr_biobert.npy", fpr)
    np.save("tpr_biobert.npy", tpr)
    with open("auc_biobert.txt", "w") as f:
        f.write(f"{roc_auc:.6f}")

    cm = confusion_matrix(true_labels, predictions)
    print(f'Confusion Matrix:\n{cm}')


if __name__ == "__main__":
    tokenizer = BertTokenizer.from_pretrained('dmis-lab/biobert-base-cased-v1.1')

    file_path = 'o8G_1.csv' 
    test_file_path = 'test_data_1.csv'  

    test_dataset = preprocess_data(test_file_path, tokenizer)
    test_dataloader = DataLoader(test_dataset, batch_size=16, shuffle=False)

    model = BertForSequenceClassification.from_pretrained('./finetuned_biobert_1')

    device = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
    model.to(device)

    evaluate_model(model, test_dataloader, device)



