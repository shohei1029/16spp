#!/bin/env zsh
set -euo pipefail

# 2016.4.11 by Shohei Nagata
# env領域のアミノ酸配列抽出，総当りblast，sim作成，mafftによるアライメントまで
# v1.0.1
#  生成されるsimファイル名のバグを修正
# 4.26
# vif 用へ。ドメイン等なしで。
# その後
# tatとか。他のタンパク質でも使えるぜ。
# --連結--
# ドメインごとの領域を対象にsimファイルの生成まで。
# ドメインはisolate HXB2を基準に。ドメインの位置情報を手動でPfamのENV HXB2のページからとってきて，toml形式とかで保存しておく。

# 5.18
# polではpiを段階的に変更してsimを構築するやつをあまり細かくやってなかったのでやる。
# polはA,B,Cの配列のみでやったblast結果がないため，全subtype対象にしてsimを構築し，そこからA,B,Cのみを切り出す.


#memo
#irisでやる設計
#stderrをログへ出力する形で実行したいね
#そのうち，swissprotのtxtファイルを読みこませるだけで，どの領域とってきて〜ってのから自動でできるようにしたいかも。でもdomain使うかchainデータ使うかとか目で見ていきたさがある。

NAME_CORE=HIV-1-gM-noRs_pol-aa
WD=`pwd`

#######################
# create sif/sim file #
#######################

mkdir -p ../analysis/sim/${NAME_CORE}

BLAST_OUT_FILE_NAME=blastp_1e-5_HIV-1-gM-noRs_pol-aa_v3.txt

NAME_CORE_ABC=HIV-1-gM-A,B,C-noRs_pol-aa

echo "creating sim files.."
#hoge.sim.txtじゃなく，hoge.txtになってしまう。。
for pi in 0 {5..6}{0,5} {7..9}{0,2,4,5,6,8}
do
#    cat ../../15spp/scripts/blast/mbs_${BLAST_OUT_FILE_NAME} | python ./v5_blast7_tosims.py -p $pi > ../analysis/sim/${NAME_CORE}/pi${pi}_${BLAST_OUT_FILE_NAME}
    perl -lane 'print "$F[0] $F[1] $F[2]" if ($F[0] =~ /.+subtype:A/ || $F[0] =~ /.+subtype:C/ || $F[0] =~ /.+subtype:B/) && ($F[1] =~ /.+subtype:A/ || $F[1] =~ /.+subtype:C/ || $F[1] =~ /.+subtype:B/)' ../analysis/sim/${NAME_CORE}/pi${pi}_${BLAST_OUT_FILE_NAME} > ../analysis/sim/${NAME_CORE_ABC}/A,B,C_pi${pi}_${BLAST_OUT_FILE_NAME}
done

jobdone create sif/sim file ${NAME_CORE} @iris, now you can visualise your network :D

# transfer files #
#nohup scp -r ../analysis/sim/${NAME_CORE} shohei@133.27.17.109:/home/shohei/bio/16spp/analysis/sim &
nohup scp -r ../analysis/sim/${NAME_CORE_ABC} shohei@133.27.17.109:/home/shohei/bio/16spp/analysis/sim &
