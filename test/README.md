VirTaK.py was executed using the following code:

```{bash, eval=FALSE, echo=TRUE}
VirTaK.py -l list.txt -d ~/db/virtak_db/VMR_MSL38_v1_complete_genomes -p ~/db/virtak_db/ -o results.txt
```

The output includes:

* results.txt
* .faa and .pfamscan files, which are further used by PanPhylo.py

Within the same working directory, PanPhylo.py was executed with the following code:

```{bash, eval=FALSE, echo=TRUE}
PanPhylo.py -l list.txt -d ~/db/virtak_db/VMR_MSL38_v1_complete_genomes -s strings.txt
```

PanPhylo directory was created
