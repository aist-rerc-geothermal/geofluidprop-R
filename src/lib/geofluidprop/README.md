地熱流体物性計算ライブラリ geofluidprop ver1.0
====

[![CMake](https://github.com/aist-rerc-geothermal/geofluidprop/actions/workflows/build_test.yml/badge.svg)](https://github.com/aist-rerc-geothermal/geofluidprop/actions/workflows/build_test.yml)

## 概要
`geofluidprop`は，主に地熱系を対象にした各種解析において必要になる流体の相状態や物性値を計算するためのC言語ライブラリです。現時点での本ライブラリの特徴は以下の通りです。
- 一般的に使用される地熱流体モデルの状態方程式を実装（現時点では，純水とH<sub>2</sub>O-NaClのみ)
- 共通インターフェイスの導入により，流体モデルの変更が用意
- C/C++等のコードから呼び出し可能（Linux, Windowsの一部環境において動作確認済み） 
- R / Pythonで使用するための簡易インターフェイスを用意 (SWIGを使用；実験的機能)
- 多倍長精度計算（一部計算式のみ）

なお，本ライブラリの開発は，研究活動の一貫として実施しており，今後も品質改良，機能拡張を行う予定です。バグの発見や機能追加の要望等がありましたら，GitHubのIssuesに投稿して頂けますと幸いです。


## コンパイル方法
本ライブラリのコンパイル済みバイナリファイルを[リリースページ](https://github.com/aist-rerc-geothermal/geofluidprop/releases)にて提供予定ですが，環境によっては動作しない，または利用可能な機能が限定される等の問題が生じる場合があります。そのため，必要に応じて，各利用者の環境で本ライブラリをソースコードからコンパイルしてください。

本ライブラリをコンパイルするためには最低限，以下のツールが必要になりますので，事前インストールをお願いします。
- C言語コンパイラ  (C11対応／GCC, Visual Studio等)
- CMake

本ライブラリのソースコードをgit cloneまたは[リリースページ](https://github.com/aist-rerc-geothermal/geofluidprop/releases)からダウンロードした後，以下のコマンドを実行することにより，ライブラリをコンパイルできます。  

Linux環境：
``` bash
cd [source-directory]
mkdir build
cd build
cmake ..
make
```

Windows環境（Visual Studio 2022を想定）：
``` bash
cd [source-directory]
mkdir build
cd build
cmake .. -G "Visual Studio 17 2022"
cmake --build . --config Release
```

コンパイルが正常に終了すると，ライブラリのバイナリファイルは `build/lib`，`build/src/Release`等に生成されます。
本ライブラリを使用してプログラム開発等を行う場合，前述のバイナリファイルに加え，`src/interface`ディレクトリに存在する以下のヘッダファイルを作業ディレクトリにコピーしてください。なお，使用する関数によってはその他のヘッダファイルもコピーする必要があります。
- `eos.h`
- `eos_type.h`
- `eos_args.h`


## 共有インターフェイスの使用例

``` C
#include <stdio.h>
#include "eos.h"

void main()
{
    EOS* eos = eos_create(EOS_TYPE_WATER_IAPWS95);

    EOS_ARGS args;
    args.p = 1e6; // [Pa]
    args.T = 273.15 + 10; // [K]
    double rho = eos_rho_pT(eos, &args); // [kg/m^3]

    printf("Water density at %g Pa and %g K is %g kg/m^3\n", args.p, args.T, rho);

    eos_free(eos);
}
```

共有インターフェイスに対応している流体モデルの一覧，及び実装している関数の詳細については， [`src/interface/eos_type.h`](src/interface/eos_type.h)，[`src/interface/eos.h`](src/interface/eos.h) をご確認ください。

また，下記のように，[`src/model`](src/model)ディレクトリ下にある低レベル関数を直接呼び出すことも可能です。
``` C
double T = 300; // [K]
double p = 0.05e6; // [Pa]
double rho = iapws95_rho_pT(p, T);
```

## SBTL機能の使用について
本ライブラリでは，一部の状態方程式について，SBTL法を用いた高速近似計算モデルを使用することが可能です。現在，IAPWS-95のみに対応しています。
当該近似計算モデルを使用するためには，[リリースページ](https://github.com/aist-rerc-geothermal/geofluidprop/releases)にある*-sbtl.zipをダウンロードする必要があります。IAPWS-95のSBTL版を使用する場合，圧縮ファイルを任意の場所に展開した後，`iapws95_sbtl_ph_table_list.txt`が含まれるディレクトリを環境変数`SBTL_IAPWS95_DIR`に設定してください。



## R，Pythonインターフェイスについて
機能は限られますが，本ライブラリをR，Python言語から呼び出すことが可能です。詳細については，[`wrapper/R/README.md`](wrapper/R/README.md)，または[`wrapper/Python/README.md`](wrapper/Python/README.md) をご確認ください。


## サードパーティライブラリの使用について
本ライブラリは以下のOSSプロジェクトの成果を使用しています。
- utest.h (https://github.com/sheredom/utest.h, Unlicense license)
- PROST (http://fluidos.etsii.upm.es/faculty/Jaime_Carpio/Fumatas_negas/PROST%20Properties%20of%20Water%20and%20Steam.htm, LGPL v2.0) (optional)
- freesteam (https://freesteam.sourceforge.net/, GPL v2) (optional)
- VTK v8.2.0 (https://github.com/Kitware/VTK, BSD 3-clause License) (optional, we use their FindR.cmake)

## 謝辞
本ライブラリが提供する機能の一部は，国立研究開発法人新エネルギー・産業技術総合開発機構（ＮＥＤＯ）の委託業務により開発されたものです。

## ライセンス
Copyright (C) 2022, National Institute of Advanced Industrial Science and Technology (AIST), All rights reserved.  

本ライブラリは，国立研究開発法人産業技術総合研究所（以下，産総研）が著作権を有しており，LGPL v3，GPL v2，または個別契約のトリプルライセンスの下で提供します。
ライセンスの詳細については，[LICENSE.md](LICENSE.md)をご確認ください。

## 免責事項
産総研は，利用者が本ライブラリ及びその派生物を使用して行う一切の行為について何ら責任を負いません。
また，本ライブラリは，予告なく，その開発，提供を中止することがあります。

## 本ライブラリを使用した成果公表について
本ライブラリを使用して得られた研究成果等を報告書，学会プロシーディング，論文等において公表される場合は，本ライブラリを使用したことについて謝辞に記載して頂けますと幸いです。

## 連絡先
本ライブラリについてご質問等がある場合は，GitHubのIssuesに投稿するか，または以下の宛先にご連絡ください。  
norihiro dot watanabe at aist.go.jp
