#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import torch
import torch.nn as nn
import torch.optim as optim
from transformers import BertTokenizer
from torch.utils.data import Dataset, DataLoader
import pandas as pd
import json
import numpy as np
from tqdm import tqdm
from model import BERT2OME  

class RNADataset(Dataset):
    def __init__(self, jsonl_path, tokenizer, max_len=43):
        self.data = []
        with open(jsonl_path, 'r') as f:
            for line in f:
                self.data.append(json.loads(line))
        self.tokenizer = tokenizer
        self.max_len = max_len

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        item = self.data[idx]
        encodings = self.tokenizer(
            item["tokenized_seq"],
            padding='max_length',
            truncation=True,
            max_length=self.max_len,
            return_tensors="pt"
        )
        chem_feat = torch.tensor(item["chem_features"], dtype=torch.float32).view(-1)
        label = torch.tensor(item["label"], dtype=torch.long)

        return {
            'input_ids': encodings["input_ids"].squeeze(0),
            'attention_mask': encodings["attention_mask"].squeeze(0),
            'chem_features': chem_feat,
            'label': label
        }

def train(model, dataloader, optimizer, criterion, device, epochs=5):
    model.train()
    for epoch in range(epochs):
        total_loss = 0
        for batch in tqdm(dataloader, desc=f"Epoch {epoch+1}"):
            input_ids = batch['input_ids'].to(device)
            attention_mask = batch['attention_mask'].to(device)
            chem_features = batch['chem_features'].to(device)
            labels = batch['label'].to(device)

            optimizer.zero_grad()
            outputs = model(input_ids=input_ids, attention_mask=attention_mask, chem_features=chem_features)
            loss = criterion(outputs, labels)
            loss.backward()
            optimizer.step()
            total_loss += loss.item()
        print(f"Epoch {epoch+1} Loss: {total_loss:.4f}")

def main():
    tokenizer = BertTokenizer.from_pretrained('bert-base-uncased')
    dataset = RNADataset("processed_data.csv", tokenizer)
    dataloader = DataLoader(dataset, batch_size=16, shuffle=True)

    device = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
    print(f"Using device: {device}")
    model = BERT2OME().to(device)

    optimizer = optim.Adam(model.parameters(), lr=1e-4)
    criterion = nn.CrossEntropyLoss()

    train(model, dataloader, optimizer, criterion, device, epochs=5)
    torch.save(model.state_dict(), "bert2ome_finetuned.pth")
    print("✅ 模型已保存为 bert2ome_finetuned.pth")

if __name__ == "__main__":
    main()


# In[ ]:


from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score

from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score, roc_curve, auc
import numpy as np

def evaluate(model, dataloader, device):
    model.eval()
    all_preds = []
    all_probs = []
    all_labels = []

    with torch.no_grad():
        for batch in dataloader:
            input_ids = batch['input_ids'].to(device)
            attention_mask = batch['attention_mask'].to(device)
            chem_features = batch['chem_features'].to(device)
            labels = batch['label'].to(device)

            outputs = model(input_ids=input_ids, attention_mask=attention_mask, chem_features=chem_features)
            probs = torch.softmax(outputs, dim=1)[:, 1]  # 概率
            preds = torch.argmax(outputs, dim=1)

            all_preds.extend(preds.cpu().numpy())
            all_probs.extend(probs.cpu().numpy())
            all_labels.extend(labels.cpu().numpy())

    acc = accuracy_score(all_labels, all_preds)
    pre = precision_score(all_labels, all_preds)
    rec = recall_score(all_labels, all_preds)
    f1 = f1_score(all_labels, all_preds)
    auc_value = roc_auc_score(all_labels, all_probs)

    print(f"✅ Evaluation results:")
    print(f"Accuracy : {acc:.4f}")
    print(f"Precision: {pre:.4f}")
    print(f"Recall   : {rec:.4f}")
    print(f"F1-score : {f1:.4f}")
    print(f"ROC AUC  : {auc_value:.4f}")

    fpr, tpr, _ = roc_curve(all_labels, all_probs)
    np.save("fpr_bert2ome.npy", fpr)
    np.save("tpr_bert2ome.npy", tpr)
    with open("auc_bert2ome.txt", "w") as f:
        f.write(f"{auc_value:.6f}")
    print("📁 Saved fpr_bert2ome.npy, tpr_bert2ome.npy, auc_bert2ome.txt")



# In[ ]:




