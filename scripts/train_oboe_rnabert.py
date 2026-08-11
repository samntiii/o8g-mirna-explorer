#!/usr/bin/env python3
"""Fine-tune multimolecule/rnabert on OBOE (Xia et al.) o8G window labels.

Figshare DOI 10.6084/m9.figshare.29634239 ships training code + CSVs but not
weights. This mirrors their rnabert/code/train.py (including one-hot labels
required by multimolecule's multilabel sequence head).

Usage (from repo root, project venv):
  python scripts/train_oboe_rnabert.py
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from sklearn.metrics import (
    accuracy_score,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
)

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA = ROOT / "third_party" / "oboe" / "rnabert" / "data"
DEFAULT_OUT = ROOT / "models" / "oboe_rnabert"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, default=DEFAULT_DATA)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--epochs", type=int, default=5)
    ap.add_argument("--batch", type=int, default=16)
    ap.add_argument("--lr", type=float, default=2e-5)
    ap.add_argument("--max-length", type=int, default=256)
    args = ap.parse_args()

    import multimolecule  # noqa: F401 — registers RNABERT
    import datasets
    from transformers import (
        AutoModelForSequenceClassification,
        AutoTokenizer,
        Trainer,
        TrainingArguments,
    )

    train_csv = args.data / "train_0.9.csv"
    valid_csv = args.data / "valid_0.9.csv"
    test_csv = args.data / "test_0.9.csv"
    if not train_csv.exists():
        raise SystemExit(f"Missing {train_csv}")

    ds = datasets.load_dataset(
        "csv",
        data_files={
            "train": str(train_csv),
            "validation": str(valid_csv),
            "test": str(test_csv),
        },
    )
    ds = ds.map(
        lambda e: {
            "sequence": [s.replace("T", "U").replace("t", "u") for s in e["sequence"]]
        },
        batched=True,
    )

    model_name = "multimolecule/rnabert"
    tokenizer = AutoTokenizer.from_pretrained(model_name, bos_token=None, eos_token=None)
    model = AutoModelForSequenceClassification.from_pretrained(model_name, num_labels=2)
    # multimolecule Criterion uses {'regression','binary','multiclass','multilabel'};
    # HF PretrainedConfig.problem_type uses different literals — only set the criterion.
    crit = model.sequence_head.criterion
    crit.problem_type = "multiclass"

    def preprocess(batch):
        enc = tokenizer(
            batch["sequence"],
            padding="max_length",
            truncation=True,
            max_length=args.max_length,
        )
        enc["labels"] = [int(float(x)) for x in batch["label"]]
        return enc

    keep = set(ds["train"].column_names)
    ds = ds.map(preprocess, batched=True, remove_columns=list(keep))
    ds.set_format(type="torch")

    def compute_metrics(pred):
        logits = pred.predictions
        if isinstance(logits, tuple):
            logits = logits[0]
        e = np.exp(logits - logits.max(axis=1, keepdims=True))
        p1 = e[:, 1] / e.sum(axis=1)
        labels = pred.label_ids
        preds = np.argmax(logits, axis=1)
        try:
            auc = float(roc_auc_score(labels, p1))
        except ValueError:
            auc = float("nan")
        return {
            "accuracy": float(accuracy_score(labels, preds)),
            "precision": float(precision_score(labels, preds, zero_division=0)),
            "recall": float(recall_score(labels, preds, zero_division=0)),
            "f1": float(f1_score(labels, preds, zero_division=0)),
            "auc": auc,
        }

    args.out.mkdir(parents=True, exist_ok=True)
    targs = TrainingArguments(
        output_dir=str(args.out / "runs"),
        num_train_epochs=args.epochs,
        per_device_train_batch_size=args.batch,
        per_device_eval_batch_size=args.batch,
        eval_strategy="epoch",
        save_strategy="epoch",
        logging_steps=50,
        learning_rate=args.lr,
        load_best_model_at_end=True,
        metric_for_best_model="auc",
        greater_is_better=True,
        report_to=[],
        save_total_limit=2,
        dataloader_pin_memory=False,
    )

    trainer = Trainer(
        model=model,
        args=targs,
        train_dataset=ds["train"],
        eval_dataset=ds["validation"],
        compute_metrics=compute_metrics,
    )
    trainer.train()
    metrics = trainer.evaluate(ds["test"])
    print("TEST metrics:", metrics)

    final = args.out / "checkpoint"
    final.mkdir(parents=True, exist_ok=True)
    trainer.model.save_pretrained(final)
    tokenizer.save_pretrained(final)
    (args.out / "test_metrics.txt").write_text(
        "\n".join(f"{k}={v}" for k, v in sorted(metrics.items())) + "\n"
    )
    print(f"Saved fine-tuned OBOE RNABERT → {final}")


if __name__ == "__main__":
    main()
