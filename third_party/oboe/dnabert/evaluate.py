#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import torch
import pandas as pd
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score, confusion_matrix
from transformers import AutoTokenizer, AutoModelForSequenceClassification
from torch.utils.data import Dataset, DataLoader

MODEL_PATH = './finetuned_dnabert_6mer'
MODEL_NAME = 'zhihan1996/DNA_bert_6'
TEST_FILE = 'test_data.csv' 
K = 6
MAX_LENGTH = 512
BATCH_SIZE = 16

def kmer_tokenize(sequence, k=6):
    sequence = sequence.upper()
    if len(sequence) < k:
        return sequence
    return ' '.join([sequence[i:i + k] for i in range(len(sequence) - k + 1)])

class TestRNADataset(Dataset):
    def __init__(self, dataframe, tokenizer, max_length=512, k=6):
        self.sequences = dataframe['sequence'].apply(lambda x: kmer_tokenize(x, k)).values
        self.labels = dataframe['label'].values
        self.tokenizer = tokenizer
        self.max_length = max_length

    def __len__(self):
        return len(self.sequences)

    def __getitem__(self, idx):
        sequence = self.sequences[idx]
        label = self.labels[idx]
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

def evaluate(model, dataloader, device):
    model.eval()
    preds = []
    true_labels = []
    probs = []

    with torch.no_grad():
        for batch in dataloader:
            input_ids = batch['input_ids'].to(device)
            attention_mask = batch['attention_mask'].to(device)
            labels = batch['labels'].to(device)

            outputs = model(input_ids, attention_mask=attention_mask)
            logits = outputs.logits
            prob = torch.softmax(logits, dim=1)[:, 1].detach().cpu().numpy()
            pred = torch.argmax(logits, dim=1)

            preds.extend(pred.cpu().numpy())
            true_labels.extend(labels.cpu().numpy())
            probs.extend(prob)

    acc = accuracy_score(true_labels, preds)
    prec = precision_score(true_labels, preds)
    rec = recall_score(true_labels, preds)
    f1 = f1_score(true_labels, preds)
    auc = roc_auc_score(true_labels, probs)
    cm = confusion_matrix(true_labels, preds)

    print(f"✅ Accuracy: {acc:.4f}")
    print(f"🎯 Precision: {prec:.4f}")
    print(f"📈 Recall: {rec:.4f}")
    print(f"📊 F1 Score: {f1:.4f}")
    print(f"🔥 AUC: {auc:.4f}")
    print("🧩 Confusion Matrix:")
    print(cm)

    from sklearn.metrics import roc_curve
    import numpy as np

    fpr, tpr, thresholds = roc_curve(true_labels, probs)
    model_tag = "DNABERT_6mer_1"  
    np.save(f"fpr_{model_tag}.npy", fpr)
    np.save(f"tpr_{model_tag}.npy", tpr)
    np.save(f"thresholds_{model_tag}.npy", thresholds)
    print(f"💾 ROC 数据保存为 fpr_{model_tag}.npy 和 tpr_{model_tag}.npy")

if __name__ == "__main__":
    print("🔍 正在加载测试数据与模型...")
    df = pd.read_csv(TEST_FILE)
    tokenizer = AutoTokenizer.from_pretrained(MODEL_NAME)
    dataset = TestRNADataset(df, tokenizer, MAX_LENGTH, K)
    dataloader = DataLoader(dataset, batch_size=BATCH_SIZE)

    device = torch.device("mps" if torch.backends.mps.is_available() else "cuda" if torch.cuda.is_available() else "cpu")
    model = AutoModelForSequenceClassification.from_pretrained(MODEL_PATH)
    model.to(device)

    evaluate(model, dataloader, device)

