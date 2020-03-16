# PME
Project Mathematical Engineering

## TODO vanaf 16/03:
#### Noodzakelijk:
- [X] done
- [ ] vullen van matrices aanpassen naar nieuwe standaard
  - [ ] main.cpp & .hpp
  - [ ] integrals.cpp & .hpp
- [ ] implementatie nonlinear solver
- [ ] alle testscenarios uittesten
  - [ ] C-waarden/plots opslaan (.m file?)
#### Minder noodzakelijk, maar zeker niet overbodig
- [ ] testrouine uitwerken (test.cpp)
  - [ ] analytisch ?
  - [ ] matlab ?
- [ ] mesh creeren dat fijner is aan de rand dan in het midden
- [ ] 

## TODO tegen verdediging
- [ ] verantwoording keuze mesh
  - [ ] grootte driehoekjes


## Useful links and references:


## Planning:
- 12/02 (1)
- 14/02 (2) ...
- 17/02 (3) Wiskundige formulering -> omzetting van probleem naar niet-lineair stelsel.
- 19/02 (4) Verificatie wiskundige formulering + Implementatie Matlab
  - Oscar:  Mesh integratie in C++
  - SibSem: Wiskundige formulering Matlab
- 21/02 (5) Implementatie C++ (componentsgewijs + componentsgewijs testen)   
  - Constructie mesh
  - Inlezen mesh
  - Aanmaken stiffness matrix
  - Non-linear solver
  - Voorstelling oplossing
  - Iteratie mesh verfijning
  (Stappen iteratief optimaliseren)
- 24/02 ( ) Implementatie C++ (componentsgewijs + componentsgewijs testen)
  Vervolg bovenstaande strategie 
- 28/02 (6) Integratie componenten
- 04/03 (7) Test op voorgestelde scenario's
- 20/03 (8) Buffer
- 25/03 (9) Evaluatie + presentatie


Huidig: 
- 28/02 wiskundige afleiding Latex -> fouten verwijderen
- 28/02 Matlab oplossing verifiëren
- 04/03 implementatie C++ afwerken (nonlin solver, data export)
        vergelijking met Matlab
- 07/03 test op voorgestelde scenario's
- 
     
        
## Notes:

## Why eigen?
- Eigen comes with built-in and fast matrix products and linear matrix decompositions, whereas Armadillo needs to link to external BLAS/Lapack libraries for that.
- Good online documentation
- Good for large square sparse linear systems https://scicomp.stackexchange.com/questions/30941/which-c-linear-algebra-library-is-probably-the-fastest-on-solving-huge-sparse

## Interessante zaken voor presentatie (aan te vullen tijdens implementatie):
