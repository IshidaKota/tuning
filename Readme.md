# ディレクトリ構造の説明
## ./run***
-fvcom　from /work/gy29/y29007/Github/fvcom442/FVCOM_source_new/fvcom(実行形式だが、当該ディレクトリで再現をするときはコピーしてあるmake.inc,makefileを用いて再度コンパイルする。)
-make.inc のコピー
-makefile
-contents.txt(inputの内容と動機を書いたファイル)

-input(ディレクトリ)
-nml
-mpi.sh
-nc(result file)

-others (eg. work.sh mpiのログ)

## ./tune
チューニング結果の可視化pythonスクリプトが入っている。
-tune.py 一連の可視化
-plot_contour.py x軸を時間,y軸を水深とした水温または塩分のコンターマップの作成
-observe_data 観測データの実測水深によるデータをシグマ層ごとに内挿しpandas.DataFrameでまとめるスクリプト
-plot_select.py
### ./tune/doc
できたfigureをtexでまとめるスクリプトたち

# 実行
runディレクトリ内のinput,nml,mpi.shを用いて計算する。
# ディレクトリ構造
run(fvcom441を用いたケース)
-fvcom from /work/gy29/y29007/Github/fvcom442/FVCOM_source/fvcom

run442(net_heat_fluxを用いた計算を目的としたケース)
-fvcom from /work/gy29/y29007/Github/fvcom442/FVCOM_source_new/fvcom

run27(FLAG_27 = -DHEATING_CALCULATEDをonにして、bulkによる熱収支の計算を目的としたもの)
-fvcom from /work/gy29/y29007/Github/fvcom442/FVCOM_source_new/fvcom
