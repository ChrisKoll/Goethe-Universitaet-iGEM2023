## localalign
Programm getestet mit Python 3.11
Ausführen über konsole mit
localalign.py promoter.fasta binding_site.fasta

Reihenfolge der Argumente ist wichtig und die erste Datei sollte die Größere sein

Standard Alignment Einstellungen sind match +5, mismacht -3, gap -3, 100 alignments ausgegeben
Standard Cluster Einstellung ist Toleranz von 20 bp auf jeder Seite

## binding
Programm getestet mit Python 3.11
Ausführen über konsole mit
binding.py promoter.fasta

Argument ist nur die Sequenz, in der gesucht werden soll alle anderen parameter können in der main Funktion geändert werden