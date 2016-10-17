
import sys

# 2016.10.17 by Shohei Nagata
# fastaファイルを読み込ませ，headerを取ってきて，そこのIDの多様性を評価する(IDの冒頭4文字とってきたときに重複なしで何個になるか).
# IDが連番とかだったら同じプロジェクトで読まれたもので非常に近い配列の可能性がある(同一個体内とか)から，そういうのを除きたいときとかに使う.

def extract_and_count_idhead(file_txt):
    idhead_s = set()
    for line in file_txt:
        line = line.rstrip()
        if '>' in line:
            idhead = line[1:5]
            idhead_s.add(idhead)
    return idhead_s

if __name__ == '__main__':
    #引数で与える形式と標準入力より与える形式の両方に対応
    if len(sys.argv) == 2:
        with open(sys.argv[1], 'r') as fh:
            idhead_s = extract_and_count_idhead(fh)
    else:
        idhead_s = extract_and_count_idhead(sys.stdin)

    print(len(idhead_s), end='')
#    print(sorted(idhead_s))
