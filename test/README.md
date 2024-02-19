VirTaK was executed using the following code:

```{bash, eval=FALSE, echo=TRUE}
VirTaK.py -l list.txt -d ~/db/virtak_db/VMR_MSL38_v1_complete_genomes -p ~/db/virtak_db/ -o virtak_results.txt
```

The output includes:
*

Within the same working directory, PanPhylo was executed with the following code:

```{bash, eval=FALSE, echo=TRUE}
PanPhylo.py -l list.txt -d ~/db/virtak_db/VMR_MSL38_v1_complete_genomes -s strings.txt
```

PanPhylo directory was created
