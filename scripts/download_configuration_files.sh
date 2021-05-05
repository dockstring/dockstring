# ExCAPE genes are the keys, and DUD-E genes are the values
declare -A excape2dude
excape2dude['ABL1']='abl1'
excape2dude['ADAM17']='ada17'
excape2dude['ADRB1']='adrb1'
excape2dude['ADRB2']='adrb2'
excape2dude['AKT2']='akt2'
excape2dude['MAOB']='aofb'
excape2dude['CASP3']='casp3'
excape2dude['DHFR']='dyr'
excape2dude['ESR2']='esr2'
excape2dude['PTK2']='fak1'
excape2dude['FGFR1']='fgfr1'
excape2dude['HMGCR']='hmdh'
excape2dude['HSP90AA1']='hs90a'
excape2dude['KIT']='kit'
excape2dude['MAPKAPK2']='mapk2'
excape2dude['MAP2K1']='mp2k1'
excape2dude['NOS1']='nos1'
excape2dude['PARP1']='parp1'
excape2dude['PDE5A']='pde5a'
excape2dude['PPARD']='ppard'
excape2dude['PGR']='prgr'
excape2dude['PTPN1']='ptn1'
excape2dude['ROCK1']='rock1'
excape2dude['AKT1']='akt1'
excape2dude['AR']='andr'
excape2dude['CDK2']='cdk2'
excape2dude['CSF1R']='csf1r'
excape2dude['ESR1']='esr1'
excape2dude['NR3C1']='gcr'
excape2dude['IGF1R']='igf1r'
excape2dude['JAK2']='jak2'
excape2dude['LCK']='lck'
excape2dude['MET']='met'
excape2dude['MMP13']='mmp13'
excape2dude['PTGS2']='pgh2'
excape2dude['PPARA']='ppara'
excape2dude['PPARG']='pparg'
excape2dude['REN']='reni'
excape2dude['ADORA2A']='aa2ar'
excape2dude['ACHE']='aces'
excape2dude['BACE1']='bace1'
excape2dude['CA2']='cah2'
excape2dude['CYP2C9']='cp2c9'
excape2dude['CYP3A4']='cp3a4'
excape2dude['HSD11B1']='dhi1'
excape2dude['DPP4']='dpp4'
excape2dude['DRD3']='drd3'
excape2dude['EGFR']='egfr'
excape2dude['F10']='fa10'
excape2dude['GBA']='glcm'
excape2dude['MAPK1']='mk01'
excape2dude['MAPK14']='mk14'
excape2dude['PLK1']='plk1'
excape2dude['SRC']='src'
excape2dude['THRB']='thb'
excape2dude['F2']='thrb'
excape2dude['KDR']='vgfr2'


# Download pdbs into 'receptors' directory
cd ~/repos/dockgym/receptors/

# Use ExCAPE name. Maybe in the supplementary info you could put a table with the correspondence
# between ExCAPE name and DUD-E name
for excape_gene in "${!excape2dude[@]}"; do
    # Echo information
    dude_gene="${excape2dude[$excape_gene]}"
    echo "ExCAPE gene: $excape_gene"
    echo "DUD-E gene: $dude_gene"
    # Download configuration file
    scp hpc_cam_cpu:/home/mg770/repos/docking/tasks/docking_benchmark/OUTPUT/create_configuration_files/"$excape_gene"_conf.txt .
done




