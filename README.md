# Stairway-To-Rainbow

Dependencies:

	- `libmpich-dev`

To execute the code run:

```
$ make
$ python3 stairwayToRainbow.py
```


As it is the code generates stepped rainbow tables with following parameters :

	- N = 2**24
	- t = 1000
	- steps positions = [750,800,850,900,950]
	- alpha = 0.95333
	- l = 2
	- nhashingnodes = 5

 
and performs :

	- nbAttack = 1000

These parameters allow to test the implementation quickly on a personal laptop.

Tables are saved in datas/ repository, informations about the precomputation phase are saved in the reccord.txt file.

Informations about the attack phase (coverage, time of the attack) are savec in output/ repository.


