# AutomaticReactionDiscovery
Using quantum chemical computation to find important reactions without
requiring human intuition.
This method help us to systematically search Reaction Mechanism Network.

## Usage
When you have configured the settings, you can use the following commands to run the program.
But this command only generate the `reactant.xyz`, `product.xyz` and `driving coordinate`.
For kinetic verification you should use the crontab and `Q-Chem` or `ORCA` with `Single-ended String Methid(SSM)` and `Pysisyphus(for irc)`.
`python ard.py input.txt reactant.xyz`
More complicate command example:
`python ard.py input.txt reactant.xyz -bonds <bonds.txt> -fixed_atoms <fixed_atoms.txt> -generations <1>`

## Preparation

 * `input.txt` : The settings
 * `reactant.xyz` : The reactant xyz which constrain all molecules
 * `bonds.txt`: For manual bond case is to specify the bonds. For cluster case is to specify missing bond in active site
 * `fixed_atoms.txt` : For catalyst system, change this may help to specify `bond_can_form`, `bond_can_break` and arrange. For openbabel and mopac constrained optimization. (xtb should also set up this because the openbabel is used before arrange)
 * `Config/` : See the following introduction

 `I will put more detailed guidelines on the wiki.`
## Config

* `orca_freq_opt_freq.lot` = Use for optimization of the irc geometry after downhill or optimization of the reactant.
Why optimize reactant? Because ARD code will have a small change of reactant during arranging the reactant/product while generating product geometry.
Calculate `Hessian` in the begining or not can be tune in here.
PAL if not necessary to be setted. Change nprocs in the database/launcher.py
Rijcosx, RIKJ or def2/J can speed up the calculation but it's a black box procedure, its use the approximation or auxiliary basis to do this.
BTW, because of `rijcosx` is a approximation method so during optimization it may cause the `imaginary(negative) frequency`

`Constraint optimization` should be setup here. Maybe like the following format index start from 0.
```
%geom
Constraints
	{C 24 C} 
	{C 25 C}
	{C 26 C}
	{C 27 C}
	end
end
```
`Partial hessian or hybrid hessian` may be setup here
```
%freq
PARTIAL_Hess 
	{24 25 26 27 28 29 30 31 32 33 34 35}
	end
end
```

* `orca_qstart` = level of theory during ssm, recommend use the `XTB2` which is the xtb `GFN2-xtb` a semi-empirical method.
* `pysisyphus_irc.yaml` = irc
Note that frozen atom in irc is achieved by isotope the Hcap atom.
But this feature is in the branch of the pysisyphus. (metadyn)
When you use the tool, please check whether this feature had been merge into master branch.
index start from 0.
```
geom:
 type: cart
 isotopes: [[24, 1000000000000000.0], [25, 1000000000000000.0], [26, 1000000000000000.0], [27, 1000000000000000.0], [28, 1000000000000000.0], [29, 1000000000000000.0], [30, 1000000000000000.0], [31, 1000000000000000.0], [32, 1000000000000000.0], [33, 1000000000000000.0], [34, 1000000000000000.0], [35, 1000000000000000.0]]
 fn: ts_geo.xyz
```
BTW, initial step of irc is from imaginary frequency so partial hessian is necessary.
But ORCA did not support all of the analytical vibrational analysis so numerical may help. Pysisyphus default is use the analytical.
```
calc:
 type: orca
 keywords: B3LYP D3BJ def2-SVP 
 blocks: "%freq PARTIAL_Hess {24 25 26 27 28 29 30 31 32 33 34 35} end end"
 pal: 8
 charge: 0
 mult: 1
 numfreq: True
```
* `frozen.txt` = frozen the atom during ssm, index start from 0.
* `qchem` = the qchem input, similiar with orca but index start from 1.
The pysisyphus did not support qchem but we usually calculate ts and irc at the same level of theory. So use qchem carefully.
Though the qchem b3lyp/def2-SVP with D3BJ correction is about the same with orca by3lyp def2-SVP D3BJ. 
(Minor energy different and gradient different. I do not compare the optimization geometry because maybe cause by different optimizer.)
* `xtb_constraint.inp` = xtb constrained optimization input. index start from 1
This is use to generate a `not bad` geometry after arranging during ard.

## Third party software modification

`Pysisyphus`
I modify the Hcap it will help frozen the Hcap during optimization while irc.
Disable the fragment during endopt is necessary because Hcap's index may have problem.
Otherwise just stop `endopt` and run a constrained optimization with other dft software.
* `pysisyphus/intcoords/update.py` see the fork in my repository.
```
        cart_step = Bt_inv_prim.T.dot(remaining_int_step)
        if nHcap != 0:
            cart_step[-(nHcap * 3):] = 0
```
`Pygsm`
version:`2745fb2417f21466ba3d33d20c6c0df11ce47fa3`
I modified the orca.py and had been pull request to the original repository. But I lost a `r` in my pull request.....
* `pygsm/level_of_theories/orca.py`
```
orcascr = 'temporcarun'
unscr ='/tmp/'+pbsID+'/'+orcasc   <<-- this should be unscr ='/tmp/'+pbsID+'/'+orcascr
```
