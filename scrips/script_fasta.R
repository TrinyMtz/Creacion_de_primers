####      Código pero que lo convierta a FASTA

library(Biostrings)

secun_prueba1 <- readDNAStringSet("extras/gallus.fasta")
secun_prueba1

### Evaluación  
pre_fw <- function (secun_prueba) {
  tipo_sec <- class(secun_prueba)
  tipo_sec[1]
  
  # Para evaluar clase de la secuencia
  if (tipo_sec[1] == "DNAStringSet") {
    # Comprobar que no es degenerado
    no_degenerado <- alphabetFrequency(secun_prueba, baseOnly=TRUE)
    
    if (no_degenerado[5] == 0) {
      # Verificar el largo de la secuencia
      longt <- width(secun_prueba)
      if ( longt < 20000) {
        # Codon de inicio: TAC -> ATG
        codon_in <- vmatchPattern ("ATG", secun_prueba)
        codon_in[[1]][1]
      } else { print("La capacidad maxima es de 20,000 nucleotidos")}}} else { print("Cambiar a DNA")}
}

pre_fw (secun_prueba1)    

######################

### FORWARD 

primer_fw <-function (inicio_codon, secun_prueba) {
  primer <- list() #Para que aparezcan opor separado del mensaje y no salgan como NULL
  detener <- inicio_codon - 20
  inicio <- 0
  ultima <- 20
  while (inicio <= detener & ultima < inicio_codon) {
    inicio <- inicio + 1  
    primer_fd <- subseq (secun_prueba, start=inicio, end=ultima) 
    ultima <- ultima + 1
    
    ##### Eveluar condiciones de primer
    # Patrones
    tripletes <- trinucleotideFrequency(primer_fd)
    patrones_malos <- tripletes [c(22, 43, 44, 16, 41, 49, 61)]
    no_hay <- c(0,0,0,0,0,0,0)
    comparacion <- all(patrones_malos == no_hay)
    
    if (comparacion == TRUE) { #Seguir evaluando el primer 
      
      # Porcentage de CG
      longt <- width(primer_fd)
      cont_cg <- letterFrequency(primer_fd, "CG")
      porc_cg <- (cont_cg / longt) * 100
      
      if (porc_cg < 60 & porc_cg > 49) {#Seguir evaluando
        #Temperatura: 55 - 65 °C
        cont_cg <- letterFrequency(primer_fd, "CG")
        cont_at <- letterFrequency(primer_fd, "AT")
        temperatura <- (4*cont_cg) + (2*cont_at)
        
        if (temperatura >54 & temperatura <66 ) {
          mensaje <- paste("Porcentaje de CG:", porc_cg, "Tm:", temperatura) 
          primer <- (list (list(primer = primer_fd, mensaje = mensaje)))
        }
      }
    }
  }
} 

primer_fw (32, secun_prueba1) #poner el valor de  interger, donde está el codón ATG    

####Guardar en un objeto 
forward_primer <- primer_fw (32, secun_prueba1)
#Obtener secuencia individual
fw_secuencia <- forward_primer[[1]]$primer
print (fw_secuencia) 
class(fw_secuencia) 


###############

### PRIMER REVERSOS
###Codones de terminación: TGA,TAG,TAA

#####Primer con TGA
revertida <- reverse(secun_prueba1)
revertida

vmatchPattern ("AGT",revertida) -> tga
tga

primer_rev_ct1 <- function (tga1,secrev) {
  primer <- list ()
  primer_encontrado <- FALSE
  detener <- tga1 - 20
  inicio <- 0
  ultima <- 20
  while (inicio <= detener & ultima < tga1 ) {
    inicio <- inicio + 1  
    primer_rv<- subseq(revertida, start=inicio, end=ultima) 
    ultima <- ultima + 1
    final<-complement(primer_rv) #la secuencia complentaria
    
    tripletes <- trinucleotideFrequency (final)#revisar patrones indeseados
    patrones_malos <- tripletes[c(22, 43, 44, 16, 41, 49, 61)]
    no_hay <- c(0,0,0,0,0,0,0)
    comparacion <- all(patrones_malos == no_hay)
    tolerancia<- sum(patrones_malos!=0)
    
    if (comparacion == TRUE | tolerancia == 1) { #Si no están esos patrones que siga evaluando
      longt <- width (primer_rv)#revisar %gc
      cont_cg <- letterFrequency(primer_rv, "CG")
      porc_cg <- (cont_cg / longt) * 100
      if (porc_cg < 60 & porc_cg > 49){ 
        cont_cg <- letterFrequency(primer_rv, "CG")
        cont_at <- letterFrequency(primer_rv, "AT")
        
        temperatura <- (4*cont_cg) + (2*cont_at) ##Tm 
        if (temperatura >54 & temperatura <66 ) {
          mensaje <- paste("Porcentaje de CG:", porc_cg, "Tm:", temperatura) 
          primer <- (list (list (primer = primer_rv, mensaje = mensaje)))
          primer_encontrado <- TRUE
        }
      }
    }
  }
  if (length(primer) == 0) {
    primer <- NULL } 
  return (primer)
} 


primer_rev_ct1 (114,revertida) #primers con TGA 

#Guardar en un objeto 
tga_rv_primer <- primer_rev_ct1(114,revertida)
#Obtener secuencia individual
tga_secuencia <- tga_rv_primer[[1]]$primer
print (tga_secuencia) 
class(tga_secuencia) #queda como un objeto en Biostrings


#####Primer con TAG

vmatchPattern ("GAT", revertida)-> tag
print (tag)

primer_rev_ct2 <- function(tag1,secrev) {
  primer <- list ()
  primer_encontrado <- FALSE
  detener <- tag1 - 17
  inicio <- 0
  ultima <- 17
  while (inicio <= detener & ultima < tag1 ) {
    inicio <- inicio + 1  
    primer_rv<- subseq(revertida, start=inicio, end=ultima) 
    ultima <- ultima + 1
    final<-complement(primer_rv)##secuencia complementaria
    
    tripletes <- trinucleotideFrequency(primer_rv)#revisar patrones indeseados
    patrones_malos <- tripletes [c(22, 43, 44, 16, 41, 49, 61)]
    no_hay <- c(0,0,0,0,0,0,0)
    comparacion <- all(patrones_malos == no_hay)
    tolerancia<-sum(patrones_malos!=0)#Más complejo que se cumplan requerimientos
    
    if (comparacion == TRUE | tolerancia==1) { #Seguir evaluando el primer 
      # Poct de CG
      longt <- width(primer_rv)
      cont_cg <- letterFrequency(primer_rv, "CG")
      porc_cg <- (cont_cg / longt) * 100
      if (porc_cg < 60 & porc_cg > 49) {##%gc
        cont_cg <- letterFrequency(primer_rv, "CG")
        cont_at <- letterFrequency(primer_rv, "AT")
        temperatura <- (4*cont_cg) + (2*cont_at) ##Tm 
        if (temperatura >54 & temperatura <66 ) {
          mensaje <- paste ("Porcentaje de CG:", porc_cg, "Tm:", temperatura) 
          primer <- ( list (list (primer = primer_rv, mensaje = mensaje)))
          primer_encontrado <- TRUE
        }
      }
    }
  } 
  if (length(primer) == 0) {
    primer <- NULL } 
  return (primer)
} 

primer_rev_ct2 (53,revertida) #sujeto a que tenga este patrón

#Guardar en un objeto 
tag_rv_primer <- primer_rev_ct2 (53,revertida)
#Obtener secuencia individual
tag_secuencia <- (tag_rv_primer[[1]]$primer)
print (tag_secuencia) 
class(tag_secuencia)



######Primers con TAA
vmatchPattern("AAT",revertida)->taa
taa

primer_rev_ct3 <- function(taa1,secrev) {
  primer <- list()
  primer_encontrado <- FALSE
  detener <- taa1 - 20
  inicio <- 0
  ultima <- 20
  while (inicio <= detener & ultima < taa1 ) {
    inicio <- inicio + 1  
    primer_rv<- subseq(revertida, start=inicio, end=ultima) 
    ultima <- ultima + 1
    tripletes <- trinucleotideFrequency(primer_rv)#revisar patrones indeseados
    patrones_malos <- tripletes[c(22, 43, 44, 16, 41, 49, 61)]
    no_hay <- c(0,0,0,0,0,0,0)
    comparacion <- all(patrones_malos == no_hay)
    if (comparacion == TRUE) { #Seguir evaluando el primer 
      # Poct de CG
      longt <- width(primer_rv)
      cont_cg <- letterFrequency(primer_rv, "CG")
      porc_cg <- (cont_cg / longt) * 100
      if (porc_cg < 60 & porc_cg > 49) {##%gc
        cont_cg <- letterFrequency(primer_rv, "CG")
        cont_at <- letterFrequency(primer_rv, "AT")
        temperatura <- (4*cont_cg) + (2*cont_at) ##Tm 
        if (temperatura >54 & temperatura <66 ) {
          mensaje <- paste("Porcentaje de CG:", porc_cg, "Tm:", temperatura) 
          primer <- (list (list (primer = primer_rv, mensaje = mensaje)))
          primer_encontrado <- TRUE
        }
      }
    }
  } 
  if (length(primer) == 0) {
    primer <- NULL } 
  return (primer)
}

primer_rev_ct3 (37,revertida)

#Guardar en un objeto 
taa_rv_primer <- primer_rev_ct3 (37,revertida)

#Obtener secuencia individual
taa_secuencia <- (taa_rv_primer[[1]]$primer)
print (taa_secuencia) 
class(taa_secuencia)


#Adición de primers válidos, que no sean NULL
primers_validos <- function (forward_primer, tga_rv_primer, tag_rv_primer, taa_rv_primer,
                              fw_secuencia, tga_secuencia, tag_secuencia, taa_secuencia) {
  lista_primers <- list()
  # Forward primer
  if (!is.null(fw_secuencia)) {
    lista_primers <- append(lista_primers, list(list(nombre = "Forward primer", secuencia = fw_secuencia)))
  }
  # TGA - Reverse primer
  if (!is.null(tga_secuencia)) {
    lista_primers <- append(lista_primers, list(list(nombre = "TGA - Reverse primer", secuencia = tga_secuencia)))
  }
  # TAG - Reverse primer
  if (!is.null(tag_secuencia)) {
    lista_primers <- append(lista_primers, list(list(nombre = "TAG - Reverse primer", secuencia = tag_secuencia)))
  }
  # TAA - Reverse primer
  if ( !is.null (taa_secuencia)) {
    lista_primers <- append (lista_primers, list (list(nombre = "TAA - Reverse primer", 
                                                       secuencia = taa_secuencia)))
    }
  return(lista_primers)
}

#Utilizar la función para solamente utilizar los primers válidos 
primers_validos (forward_primer, tga_rv_primer, tag_rv_primer, taa_rv_primer,
                 fw_secuencia, tga_secuencia, tag_secuencia, taa_secuencia)
#Guardar en un objeto
lista_primers <- primers_validos (forward_primer, tga_rv_primer, tag_rv_primer, taa_rv_primer,
                                     fw_secuencia, tga_secuencia, tag_secuencia, taa_secuencia)

#Generar listas por de secuencias y nombres para el archivo FASTA

escribir_fasta <- function (lista_primers) {
  lista_secuencias <- list ()
  lista_nombres <- list ()   
  
  for (primer in lista_primers) {
    lista_secuencias <- c(lista_secuencias, primer$secuencia)
    lista_nombres <- c(lista_nombres, primer$nombre)
  }
  
  if (length(lista_secuencias) > 0 | length(lista_nombres) > 0) {
      library(seqinr)
      write.fasta(sequences = lista_secuencias, 
                  names = lista_nombres, 
                  nbchar = 80, 
                  file.out = "scrips/resultados/primers2.fasta")
      print("Ver primers en carpeta de resultados")
    }
  }

escribir_fasta (lista_primers)












