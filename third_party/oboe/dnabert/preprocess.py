#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import torch
from torch.utils.data import Dataset, DataLoader
from transformers import AutoTokenizer

K = 6
MAX_LENGTH = 512
BATCH_SIZE = 16
FILE_PATH = 'o8G_1.csv'

def kmer_tokenize(sequence, k=6):
    sequence = sequence.upper()
    if len(sequence) < k:
        return sequence
    return ' '.join([sequence[i:i + k] for i in range(len(sequence) - k + 1)])

class RNADataset(Dataset):
    def __init__(self, dataframe, tokenizer, max_length=512, k=6):
        self.tokenizer = tokenizer
        self.max_length = max_length
        self.k = k
        self.sequences = dataframe['sequence'].apply(lambda x: kmer_tokenize(x, k)).values
        self.labels = dataframe['label'].values

    def __len__(self):
        return len(self.sequences)

    def __getitem__(self, index):
        sequence = self.sequences[index]
        label = self.labels[index]

        if index == 0:  # 只输出一次
            print("🔬 第一个样本的 k-mer 序列：", sequence)
            print("🔬 Tokenizer 输出 token:", self.tokenizer.tokenize(sequence))
            
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

def load_dataloader():
    print("🔄 正在加载数据...")
    data = pd.read_csv(FILE_PATH)
    tokenizer = AutoTokenizer.from_pretrained('zhihan1996/DNA_bert_6')
    dataset = RNADataset(data, tokenizer, MAX_LENGTH, K)
    dataloader = DataLoader(dataset, batch_size=BATCH_SIZE, shuffle=True)
    print("✅ 预处理完成，DataLoader 已就绪！")
    return dataloader

if __name__ == "__main__":
    dataloader = load_dataloader()
    for batch in dataloader:
        print({k: v.shape for k, v in batch.items()})
        break  # 只查看第一个 batch 的维度

