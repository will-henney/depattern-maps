D=../OrionMuse/LineMaps
for f in $D/mean-*-????.fits; do
    echo "Processing  $f"
    python de-pattern-extreme.py $f
done
