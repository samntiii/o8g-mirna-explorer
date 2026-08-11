#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import torch
from model import BERT2OME
from transformers import BertTokenizer
from torch.utils.data import DataLoader
from train import RNADataset, evaluate  

def main():
    tokenizer = BertTokenizer.from_pretrained('bert-base-uncased')
    dataset = RNADataset("processed_data.csv", tokenizer)
    dataloader = DataLoader(dataset, batch_size=16, shuffle=False)

    device = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
    model = BERT2OME().to(device)
    model.load_state_dict(torch.load("bert2ome_finetuned.pth", map_location=device))
    print("✅ 模型已加载，开始评估...")
    
    evaluate(model, dataloader, device)

if __name__ == "__main__":
    main()


# In[ ]:




