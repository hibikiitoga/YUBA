# YUBA
Program for generating vtk file from simulation result

##入力ファイル形式
シミュレーション結果の座標や結合情報をタグ付けして出力指定ください。
“#” から始まる行 はタグを取り扱う行とみなされます。
現在対応している型は Bead, Bond, Coordinate、Triangle, Cell です。
タグは型情報とステップ数を表す正の整数で構成されます (Cell 型を除く)。
区切り文字は半角スペースもしくはタブです。
今後、簡便のため整数値を “int”、浮動小数点数を “double” と呼びます。
(double と書かれていますが、単精度と倍精度に対応しています。)

