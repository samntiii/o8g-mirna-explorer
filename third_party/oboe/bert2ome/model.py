#!/usr/bin/env python
# coding: utf-8

# In[1]:

from transformers import BertModel
import torch
import torch.nn as nn

class BERT2OME(nn.Module):
    def __init__(self, bert_model='bert-base-uncased'):
        super(BERT2OME, self).__init__()
        self.bert = BertModel.from_pretrained(bert_model)

        self.cnn = nn.Sequential(
            nn.Conv2d(1, 32, kernel_size=(3, 3), padding=1),
            nn.ReLU(),
            nn.MaxPool2d((2, 2)),
            nn.Dropout(0.3),
            nn.Conv2d(32, 64, kernel_size=(3, 3), padding=1),
            nn.ReLU(),
            nn.MaxPool2d((2, 2)),
            nn.Flatten()
        )

        with torch.no_grad():
            dummy = torch.zeros(1, 1, 43, 768)
            cnn_out = self.cnn(dummy)
            cnn_output_size = cnn_out.view(1, -1).shape[1]

        self.fc = nn.Sequential(
            nn.Linear(cnn_output_size + 123, 256),
            nn.ReLU(),
            nn.Dropout(0.3),
            nn.Linear(256, 2)
        )

    def forward(self, input_ids, attention_mask, chem_features):
        outputs = self.bert(input_ids=input_ids, attention_mask=attention_mask)
        x = outputs.last_hidden_state.unsqueeze(1)  # shape: [B, 1, 43, 768]
        cnn_out = self.cnn(x)
        combined = torch.cat((cnn_out, chem_features), dim=1)
        return self.fc(combined)
