---
uuid: 3c6ce746-bee5-461d-afca-17969c0281fd
title: matplotlibでcolorbarのみplotし、ある値がどの色になるのか判定する
description: matplotlibを使っていて、colorbarだけ作りたいとき、そして作ったcolorbarに対して、ある値がどの色になるのか知りたい、というニッチな状況への対応策
lang: ja
category: python
created_at: 2020/09/26
updated_at: 2021/09/29
---
## TL;DR

matplotlibを使っていて、colorbarだけ作りたい！そして、何らかの値がそのcolorbarのどの色になるのか知りたい！というようなことがあります。

**ex)**
何らかのSVGがあって、それに値に応じた色をつけたい、そしてカラーバーも欲しい

 > 
 > Python 3.7.4

## やり方

[matplotlib.colorbar.ColorbarBase](https://matplotlib.org/3.3.1/api/colorbar_api.html#matplotlib.colorbar.ColorbarBase)を使います。また、カラーバーの値の範囲を決める、ある値がどの色になるかを決める際に、[matplotlib.colors.Normalize](https://matplotlib.org/3.3.1/api/_as_gen/matplotlib.colors.Normalize.html)を使います。

まず範囲を決めます。

````python
import matplotlib as mpl
import matplotlib.pyplot as plt

print(mpl.__version__)
# '3.2.2'

vmin = -10
vmax = 10

norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
````

カラーバーを書きます。[matplotlib.pyplot.get_cmap](https://matplotlib.org/3.3.1/tutorials/colors/colormaps.html)でcolormapの情報を持ってきます。範囲を決める際に、先程用意したnormを用います。

saveするときが少し注意が必要で、`bbox_inches="tight"`をオプションで指定しないとticksや、label情報が消えます。

````python
fig, ax = plt.subplots(figsize=(1, 5))
cmap = plt.get_cmap("Wistia")
cbar = mpl.colorbar.ColorbarBase(
    ax=ax,
    cmap=cmap,
    norm=norm,
    orientation="vertical",
    label="sample",
)

plt.savefig("sample_colormap.png", bbox_inches="tight")
````

<amp-img src="/images/colorbar_sample.png" height="20rem" width="8rem" alt="sample_colorbar" />

対応するrgbaカラーを取得します。

````python
norm_value = norm(5)
rgba = cmap(norm_value)
print(rgba)
# (0.9998615916955017, 0.6259284890426758, 0.0, 1.0)
````