```text 
                          _    _    ____    _____     _____   ______    _____ 
                         | |  | |  / __ \  |  __ \   / ____| |  ____|  / ____|
                         | |__| | | |  | | | |__) | | (___   | |__    | (___  
                         |  __  | | |  | | |  _  /   \___ \  |  __|    \___ \ 
                         | |  | | | |__| | | | \ \   ____) | | |____   ____) |
                         |_|  |_|  \____/  |_|  \_\ |_____/  |______| |_____/ 
                                                                              
                                                                                                                                                            
                         +m:y:-`                                       `-:y:m+                         
                 .:::odmMMMNMms:-                                   -:smMNMMMmdo:::.                   
                 `+mMMNMMMMNNNd+:                                   :+dNNNMMMMNMMm+`                   
               \shdNMNymMMMMshMmdNmd:                           :dmNdmMhsMMMMmyNMNdhs/                 
              ooosMMMNMMMNmNhNMMMNmhd-                         -dhmNMMMNhNmNMMMNMMMsooo                
             `-\+yhMMmMMMhhMhss+\+o+.                           .+o+/+sshMhhMMMmMMhy+/-`               
                 yMNMhMMMMNMs                                           sMNMMMMhMNMy                   
               :ydmMMhMMMMy+M\                                         /M+yMMMMhMMmdy:                 
              .+s+omMhMMMMMy+Ms`                                     `sM+yMMMMMhMmo+s+.                
                +ymNMNmMMMMMmdMm+`                                 `+mMdmMMMMMmNMNmy+                  
  :-           \\oymMMmMMMMMMMmNMd-       \yys:       :syy/       -dMNmMMMMMMMmMMmyo//           -:    
 - -y.          .\oymMMMMMMMMhddhNMNs..+yNMMMMm-     -mMMMMNy+..sNMNhddhMMMMMMMMmyo/.          .y- -   
 +o.-N.       `-:+yhdNMMMMMsMMMMMMMMMMmMMMmd\msd     dsm/dmMMMmMMMMMMMMMMsMMMMMNdhy+:-`       .N-.o+   
: :NsMo           oyo+-hMMMMMMMMMMMMMM+NMd+`.+mM\   /Mm+.`+dMN+MMMMMMMMMMMMMMh-+oyo           oMsN: :  
\y.yMM:          `.   :hoMMMMMMMMMMMMMMmMNMMMMmsN` `NsmMMMMNMmMMMMMMMMMMMMMMoh:   .`          :MMy.y/  
 oNMMh               oh+MMMMMMMMMMMMMMMMNdNNyNMdms smdMNyNNdNMMMMMMMMMMMMMMMM+ho               hMMNo   
 oMMM-             \dyhMMMMMMMMMMMMMMMMNs\. `mmsMm mMsmm` ./sNMMMMMMMMMMMMMMMMhyd/             -MMMo   
 dNym           -oNMMMMMMMMMMMMMMMMNoo\.    odh-mN Nm-hdo    ./ooNMMMMMMMMMMMMMMMMNo-           myNd   
`Mdsh        `+mMMMMMMMMMMMMMMMh+Mh-       -NM\ `. .` /MN-       -hM+hMMMMMMMMMMMMMMMm+`        hsdM`  
`NN:N.      oNMMMMMMNoMMMMMMMhoyh:       .oNMd`       `dMNo.       :hyohMMMMMMMoNMMMMMMNo      .N:NN`  
 \Mm\do`  .mMMMMMMMMMMMMMMMmhyo.         NNM+           +MNN         .oyhmMMMMMMMMMMMMMMMm.  `od/mM/   
  .smNMMdymMMMMMMMMNmMMMMMhy:            sy-             -ys            :yhMMMMMmNMMMMMMMMmydMMNms.    
    -NhMMdMhmMMMMMMMyMMMhhMM+                                           +MMhhMMMyMMMMMMMmhMdMMhN-      
    \yhNd-hNMMMMMMMMmymhmMMd                                             dMMmhmymMMMMMMMMNh-dNhy/      
    :-sm\-`dMMMMMMMMhh\MMMm.                                             .mMMM/hhMMMMMMMMd`-/ms-:      
      `yy  `oNMMMMMMMmsMMm.                                               .mMMsmMMMMMMMNo`  yy`        
        y`   `omMMMMMMsMy`    -\:\+++o.                       .o+++/:/-    `yMsMMMMMMmo`   `y          
                sMMhmoNNhysyhmMMmNMMho`                       `ohMMNmMMmhysyhNNomhMMs                  
                oMMN:yMmyo\-.``                                       ``.-/oymMy:NMMo                  
               oNMN.                                                             .NMNo                 
               -sNmo.             High Order Spectral Element Solver            .omNs-                 
                  \dyo.                                                       .oyd/                    
                   `oMNd-      2D Compressible Navier-Stokes equations      -dNMo`                     
                     :dMNy`                                               `yNMd:                       
                      -ymms                                               smmy-

```

***This is a two dimensional discontinuous Galerkin spectral element method solver for Navier-Stokes equations.***

This code has been developed in the Madrid Technical University.

Installation:
	
1. Replace make.inc.example to make.inc with the paths of your system.
2. Run make. Several options are available:
	* Solver: Euler (make Euler) or Navier-Stokes (make all)
	* Mode: DEBUG/RELEASE
	* Compiler: gfortran/ifort

3. Run the three benchmark tests:
	* make runcyl. Flow around a circle
	* make runvortex. Taylor vortex problem.
	* make runchan. Channel.
