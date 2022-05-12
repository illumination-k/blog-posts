---
title: 非モデル生物でOrg.db (AnnotationDbi) を自作する
description: NGS解析において、GO Enrichment解析を始めとしたEnrichment解析はデータを生物学的に解釈するために重要な手法だが、非モデル生物ではAnnotationDbiがない場合が多いので自作する方法を説明する。
---

## TL;DR

NGS解析において、GO Enrichment解析を始めとしたEnrichment解析はデータを生物学的に解釈するために重要な手法ですが、非モデル生物で行うのはAnnotationDbiなどのパッケージが充実しておらず少し敷居が高いです。そこで、[AnnotationForge]()を使ってAnnotationDbi (Org.db) を自作し、GO Enrichment解析を行います。

ちなみにAnnotationDbiの中身はSQLite3なのでパッケージに頼らなくてもどうにかなりますが、簡単なのでAnnotationForgeを利用します。

今回は、_Marchantia polyrmopha_ のデータを使用します。_Marchantia polymorpha_ は和名ではゼニゴケで、日本で活発に研究されている植物の1つです。

## 1. アノテーションファイルの準備

まず、データベースを作成するためのアノテーションを準備します。今回の目的はGO Enrichment解析なので、GO annotationが必要です。

ゼニゴケの場合は、[marpolbase](https://marchantia.info)というサイトでアノテーションが管理されているのでこちらを利用します。

もしアノテーションが利用できない場合は、InterproscanやBlast2GOなどを利用してアノテーションしてください。

gene idとgo term idが1対1対応しているようなファイルを作成します。

```csv:title=gene2go.csv
Mp1g00070,GO:0006281
Mp1g00070,GO:0006974
Mp1g00080,GO:0004853
Mp1g00080,GO:0006779
Mp1g00120,GO:0004385
Mp1g00120,GO:0006163
Mp1g00140,GO:0003723
Mp1g00180,GO:0006520
Mp1g00180,GO:0016491
Mp1g00180,GO:0016639
```

結構データが汚いので、Pythonで綺麗にします。

データを見ると、以下のような構造になっています。なので、まずは転写因子とアノテーション列を分割したあと、アノテーション部分のGO Termを持ってきます。

```
transcripts_id \t type1:description;annotation_type2;...
```

```python
gene2go = {}

with open("MpTak_v6.1_func_annotation_1line.tsv") as f:
    for line in f:
        line = line.strip()
        elems = line.split("\t")

        if len(elems) != 2:
            continue

        tx, annotations = elems
        gene = tx.split(".")[0]

        # only use a first transcript
        if gene in gene2go:
            continue

        gene2go.setdefault(gene, [])

        for annotation in annotations.split(";"):
            annotation = annotation.strip()
            if not annotation.startswith("GO:"):
                continue

            elems = annotation.split(":")
            go_term = ":".join(elems[:2])
            gene2go[gene].append(go_term)

with open("gene2go.csv", "w") as w:
    lines = []
    for gene, go_terms in gene2go.items():
        for go_term in go_terms:
            lines.append(f"{gene},{go_term}")
    lines = list(sorted(set(lines)))
    w.write("\n".join(lines))
```

## 2. 実行環境の準備

ローカルのRを使っている人は以下をコンソールで実行してください。

```R
install.packages("BiocManager")
BiocManager::install(c("AnnotationForge", "GO.db"), ask = FALSE)
```

Dockerを使いたい場合はこちらをビルドしておいてください。

```docker
FROM rocker/r-base

RUN apt-get update --no-install-recomentds --fix-missing -y && \
    apt-get install -y --no-install-reccomentds libxml2-dev libcurl4-openssl-dev openssl r-cran-openssl libssl-dev

RUN Rscript -e 'install.packages("BiocManager")' && \
    Rscript -e 'BiocManager::install(c("AnnotationForge", "GO.db"), ask = FALSE)'

WORKDIR /local
```

## 3. Rの実行

```r
library(AnnotationForge)

geneInfo <- read.csv("gene_info.csv", row.names = NULL)
gene2go <- read.table("gene2go.csv", row.names = NULL, header = TRUE, sep = ",")

makeOrgPackage(gene_info=geneInfo, go=gene2go,
               version="0.1",
               maintainer="illumination-k <illumination.k.27@gmail.com>",
               author="illumination-k <illumination.k.27@gmail.com>",
               outputDir = ".",
               tax_id="3197",
               genus="Marchantia",
               species="polymorpha",
               goTable='go')
```
