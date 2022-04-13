file=$1

sed -n '/Coverage/,$p' $file | sed 's/Coverage/#Coverage/' > $(basename $file .divsum).cov
sed '/Coverage/q' $file | grep -v "Coverage" | sed '1,6 {s/^/#/}' | sed 's/#Class/Class/' > $(basename $file .divsum).kim



