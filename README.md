# Creacion_de_primers
 
## Primers

### ¿Qué son?
Los primes o cebadores son ***secuencias cortas de nucleotidos*** que marcan el punto de partida para la ***sintesis del DNA por la DNA polimerasa***. En una PCR, permiten seleccionar la región de DNA que se busca amplificar o copiar. 

Se diseña un par de primers **_Forward (Fw)_**, complementario a la secuencia 3'-5' en el extremo 3', y **_Reverse (Rev)_**, complementario a la secuencia 5'-3' en el extremo 3'. 

### Caracteristicas
- Son de cadena sencilla.
- Generalmente son oligonucleotidos de más de 18 nucleotidos. 
- Deben ser complementarios a la secuencia blanco. 
- Son especificos para un sitio.
 
## Condiciones para su uso en PCR
1. Deben tener una longuitud de 17-24 pb.
2. Deben de tener una temperatura de fusión (Tm) de 50-65°C.
3. La Tm entre ambos primers no debe diferir en más de 5°C.
4. Que tengan un continido de G|C de 50-60%.
5. Que no presenten repeticiones de nucleotidos >4 o secuencias como:

            GGG, CCC, GGT, ATT, CGA, TAA O TTA. 

para evitar la formación de horquillas (hairpins).

6. Primers deben ser distintos entre sí, para evitar la formación de dímeros.

### Tm
Es la temperatura a la cual una parte de la secuencia de DNA se encuentra en forma 
monocatenaria, de fomrma que puede interactuar con secuencias complementarias, como 
los primers, de modo que se hibridan con el DNA. 

Se calcula a partir de la ***_Regla de Wallace_***, que se basa en el contenido de GC y AT de 
la secuencia:

                     Tm: 4GC + 2AC 

## Codigo
Se generaron una serie de funciones que permitan obtener los primers, forward y reverse, de una secuencia de DNA. Para poder generar y usar estas funciones es necesario tener instalado y cargar el paquete de Biostrings. 
El [codigo para generar los primers](scrips/evaluar.R) se encuentra disponible en este mismo repositorio.


### Limitantes
- Solo acepta secuencias de DNA menores a 20,000 nucleotidos.
- Unicamente util para secuencias que tienen el codon de inicio TAC.
- Si la secuencia disponible para hacer el primer es muy corta y/o presenta muchas secuencias que propician los hairpins, no los genera. 
