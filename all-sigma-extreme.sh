D=../OrionMuse/LineMaps
for f in $D/sigma-*-????.fits; do
    echo "Processing  $f"
    python de-pattern-extreme.py $f
done
