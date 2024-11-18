#   PRIMERS

library(Biostrings)

secun_prueba <- readDNAStringSet("extras/sequence.fasta")
secun_prueba

tipo_sec <- class(secun_prueba)
tipo_sec[1]

# Para evaluar tp de secuencia
if (tipo_sec[1] == "DNAStringSet") {
  print("correcto")
} else {
  print("Convierte tu secuencia a DNA")
}

# Comprobar que no es degenerado
no_degenerado <- alphabetFrequency(secun_prueba, baseOnly=TRUE)
no_degenerado[5]

if (no_degenerado[5] == 0) {
  print("correcto")
} else {
  print("Convierte tu secuencia a DNA")
}

# Para encontar la secuencia anterior al codon de inicio
# Verificar el largo de la secuencia
longt <- width(secun_prueba)

if ( longt < 20000) {
  print("correcto")
} else {
  print("muy larga")
}

# Codon de inicio: TAC -> ATG
codon_in <- vmatchPattern("ATG", secun_prueba)
print(codon_in[[1]][1]) 

# Antes de la sec. inicio
### secun_forward <- subseq(secun_prueba, start=1, end=92)
### secun_forward

42

## Generador de secuencias. Desplazamiento del marco de lect de n+1 a 20
detener <- 
inicio <- 0
ultima <- 20
while (inicio <= detener) {
  inicio <- inicio + 1

  primer_fd <- subseq(secun_prueba, start=inicio, end=ultima) 
  
  print(primer_fd)
  
  ultima <- ultima + 1
  
  ##### Eveluar condiciones de primer 
  
}

## Generador de secuencias. Desplazamiento del marco de lect de n+1 a 18
detener <- 74
inicio <- 0
ultima <- 18
while (inicio <= detener & ultima < 92 ) {
  inicio <- inicio + 1
  
  primer_fd <- subseq(secun_prueba, start=inicio, end=ultima) 
  
  print(primer_fd)
  
  ultima <- ultima + 1
  
  ##### Eveluar condiciones de primer 
  
}

## Generador de secuencias. Desplazamiento del marco de lect de n+1 a 24
detener <- 68
inicio <- 0
ultima <- 24
while (inicio <= detener) {
  inicio <- inicio + 1
  
  primer_fd <- subseq(secun_prueba, start=inicio, end=ultima) 
  
  print(primer_fd)
  
  ultima <- ultima + 1
  
  ##### Eveluar condiciones de primer 
}

#### REQUISISTOS

# Porcentaje de CG


# Patrones q no debe haber





