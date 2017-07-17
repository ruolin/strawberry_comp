awk '{print $1}' data/tx_exp.mat > data/tname.lst
./gff-select.py ~/Research/Annotations/hsapiens/gencode.v19.annotation.gff3 -l data/tname.lst > truth.gff
