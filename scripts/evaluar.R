###################
library(Biostrings)
install.packages ("seqinr)
                  
# Lectura de la secuencia
secun_prueba1 <- readDNAStringSet("Creacion_de_primers/extras/gallus.fasta")
secun_prueba1
                   ############################
####################     PRIMERS FORWARD      ############################
                   ############################


#   IDENTIFICACION DE LA SEC CODIF (CDS)

pre_fw <- function(secun_prueba) {
  tipo_sec <- class(secun_prueba)
  tipo_sec[1]
  
  # Para evaluar tipo de secuencia
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
        print(codon_in)
      } else { print("La capacidad maxima es de 20,000 nucleotidos")}}} else { print("Cambiar a DNA")}
  }
    

pre_fw(secun_prueba1) #Si el patrón está muy cerca, porque está más de una vez, intenete con otras posiciones

####################

#  GENERADOR DE PRIMERS FORWARD

primer_fw <-function (inicio_codon, secun_prueba, ultima) {
  fw_primers <- list()
  detener <- inicio_codon - ultima
  inicio <- 0
  while (inicio <= detener & ultima < inicio_codon ) {
          inicio <- inicio + 1  
          primer_fd <- subseq (secun_prueba, start=inicio, end=ultima) 
          ultima <- ultima + 1
          
          final<-complement(primer_fd)
          ##### Eveluar condiciones de primer
          # Patrones

          tripletes <- trinucleotideFrequency(final)
          patrones_malos <- tripletes[c(22, 43, 44, 16, 41, 49, 61)]
          no_hay <- c(0,0,0,0,0,0,0)
          comparacion <- all(patrones_malos == no_hay)
          
          if (comparacion == TRUE) { #Seguir evaluando el primer 

            # Poct de CG: 50-60 %
            longt <- width(final)
            cont_cg <- letterFrequency(final, "CG")
            porc_cg <- (cont_cg / longt) * 100
            
            if (porc_cg < 61 & porc_cg > 49) {#Seguir evaluando
              #Temperatura: 55 - 65 °C
              cont_cg <- letterFrequency(final, "CG")
              cont_at <- letterFrequency(final, "AT")
              temperatura <- (4*cont_cg) + (2*cont_at)
              
              if (temperatura >54 & temperatura <66 ) {
                fw_primers <- append (fw_primers, list (list (primer = final, porc_cg = porc_cg, temperatura = temperatura)))
              }
              }
          }
          } 
  if (length(fw_primers) > 0) {
    for (primer in fw_primers) {
      print("Forward")
      print(primer$primer)
      print(paste("Porcentaje de CG: ", primer$porc_cg))
      print(paste("Tm: ", primer$temperatura, "°C"))
    }
  } else { 
    fw_primers <- NULL
  }
  return (fw_primers)
  } 



primer_fw (32, secun_prueba1, 18) # Colocar el valor del start , secuencia, long del primer. 


                   ############################
####################     PRIMERS REVERSE      ############################
                   ############################

### Codones de terminación: TGA,TAG,TAA

####################
# Primer con TGA
revertida <-reverse(secun_prueba1)
revertida

vmatchPattern ("AGT",revertida) -> tga 
tga

primer_rev_ct1<-function(tga1,secrev,ultima) {
  tga_rv_primers <- list ()
  detener <- tga1 - ultima

  inicio <- 0
  while (inicio <= detener & ultima < tga1 ) {
    inicio <- inicio + 1  
    primer_rv<- subseq(revertida, start=inicio, end=ultima) 
    ultima <- ultima + 1
    final<-complement(primer_rv) #la secuencia complentaria
    
    tripletes <- trinucleotideFrequency(final)#patrones que favorecen horquillas,
    #dímeros..
    patrones_malos <- tripletes[c(22, 43, 44, 16, 41, 49, 61)]
    no_hay <- c(0,0,0,0,0,0,0)
    comparacion <- all(patrones_malos == no_hay)

    if (comparacion == TRUE) { 
      longt <- width(primer_rv)#revisar %gc
      cont_cg <- letterFrequency(primer_rv, "CG")
      porc_cg <- (cont_cg / longt) * 100
      
      if (porc_cg < 61 & porc_cg > 49){ 
        cont_cg <- letterFrequency(primer_rv, "CG")
        cont_at <- letterFrequency(primer_rv, "AT")
        temperatura <- (4*cont_cg) + (2*cont_at) ##Tm 
        
        if (temperatura >54 & temperatura <66 ) {
          tga_rv_primers <- append (tga_rv_primers, list (list (primer = primer_rv, porc_cg = porc_cg, temperatura = temperatura)))
        }
      }
    }
  }
  if (length(tga_rv_primers) > 0) {
    for (primer in tga_rv_primers) {
      print("reverse tga")
      print(primer$primer)
      print(paste("Porcentaje de CG: ", primer$porc_cg))
      print(paste("Tm: ", primer$temperatura, "°C"))
    }
  } else { 
    tga_rv_primers <- NULL
  }
  return (tga_rv_primers)
}  

primer_rev_ct1(114,revertida, 20) # Primers con TGA, puedes, jugar con los valores de 
# Start del objeto TGA, considerando tu región codificante (CDS)



####################
# Primer con TAG

vmatchPattern ("GAT",revertida)->tag
tag


primer_rev_ct2<-function(tag1,secrev,ultima) {
  tag_rv_primers <- list ()
  detener <- tag1 - ultima
  inicio <- 0
  while (inicio <= detener & ultima < tag1 ) {
    inicio <- inicio + 1  
    primer_rv<- subseq(revertida, start=inicio, end=ultima) 
    ultima <- ultima + 1
    final<-complement(primer_rv)##secuencia complementaria
    
    tripletes <- trinucleotideFrequency(primer_rv)#revisar patrones indeseados
    patrones_malos <- tripletes [c(22, 43, 44, 16, 41, 49, 61)]
    no_hay <- c(0,0,0,0,0,0,0)
    comparacion <- all(patrones_malos == no_hay)
    
    if (comparacion == TRUE) { 
      
      longt <- width(primer_rv)#%gc
      cont_cg <- letterFrequency(primer_rv, "CG")
      porc_cg <- (cont_cg / longt) * 100
      
      if (porc_cg < 61 & porc_cg > 49) {##%gc
        cont_cg <- letterFrequency(primer_rv, "CG")
        cont_at <- letterFrequency(primer_rv, "AT")
        temperatura <- (4*cont_cg) + (2*cont_at) ##Tm 
        
        if (temperatura >54 & temperatura <66 ) {
          tag_rv_primers <- append (tag_rv_primers, list(list (primer = primer_rv, porc_cg = porc_cg, temperatura = temperatura)))
        }
      }
    }
  }
  if (length(tag_rv_primers) > 0) {
    for (primer in tag_rv_primers) {
      print("reverse tag")
      print(primer$primer)
      print(paste("Porcentaje de CG: ", primer$porc_cg))
      print(paste("Tm: ", primer$temperatura, "°C"))
    }
  } else { 
    tag_rv_primers <- NULL
  }
  return (tag_rv_primers)
} 


primer_rev_ct2(222,revertida, 20)#reverse TAG,sujeto a que tenga este patrón,puedes
#jugar con los valores de start del objeto tga, considerando tu región codificante (CDS)


####################
# Primers con TAA
vmatchPattern("AAT",revertida)->taa
taa

primer_rev_ct3<-function(taa1,secrev,ultima) {
  taa_rv_primers <- list ()
  detener <- taa1 - ultima
  inicio <- 0
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
      
      if (porc_cg < 61 & porc_cg > 49) {##%gc
        cont_cg <- letterFrequency(primer_rv, "CG")
        cont_at <- letterFrequency(primer_rv, "AT")
        temperatura <- (4*cont_cg) + (2*cont_at) ##Tm 
        
        if (temperatura >54 & temperatura <66 ) {
          taa_rv_primers <- append (taa_rv_primers, list(list (primer = primer_rv, porc_cg = porc_cg, temperatura = temperatura)))
        }
      }
    }
  } 
  if (length(taa_rv_primers) > 0) {
    for (primer in taa_rv_primers) {
      print("reverse taa")
      print(primer$primer)
      print(paste("Porcentaje de CG: ", primer$porc_cg))
      print(paste("Tm: ", primer$temperatura, "°C"))
    }
  } else { 
    taa_rv_primers <- NULL
  }
  return (taa_rv_primers)
}

primer_rev_ct3(37,revertida, 20) #reverse TAG,sujeto a que tenga este patrón,puedes
#jugar con los valores de start del objeto tga, considerando tu región codificante (CDS)


#Objetos 
primer_fw (32, secun_prueba1, 18) -> fw_primers
print(fw_primers)
length(fw_primers)

primer_rev_ct1 (114, revertida, 20) -> tga_primers
print (tga_primers)
length(tga_primers)

primer_rev_ct2 (222, revertida, 20) -> tag_primers
print(tag_primers)
length(tag_primers)

primer_rev_ct3 (37, revertida, 20) -> taa_primers
print(taa_primers)
length(taa_primers)

#Adición de primers válidos, que no sean NULL
primers_validos <- function(fw_primers, tga_primers, tag_primers, taa_primers) {
  lista_primers <- list()
  #forward 
  if (!is.null(fw_primers)) {
    for (primer in fw_primers) {
      lista_primers <- c(lista_primers, list(list(nombre = "Forward primer", secuencias = primer$primer)))
    }
  }
  #TGA reverso
  if (!is.null(tga_primers)) {
    for (primer in tga_primers) {
      lista_primers <- c(lista_primers, list(list(nombre = "TGA - Reverse primer", secuencias = primer$primer)))
    }
  }
  #TAG reverso
  if (!is.null(tag_primers)) {
    for (primer in tag_primers) {
      lista_primers <- c(lista_primers, list(list(nombre = "TAG - Reverse primer", secuencias = primer$primer)))
    }
  }
  #TAA reverso
  if (!is.null(taa_primers)) {
    for (primer in taa_primers) {
      lista_primers <- c(lista_primers, list(list(nombre = "TAA - Reverse primer", secuencias = primer$primer)))
    }
  }
  return(lista_primers)
}


#Utilizar la función para solamente utilizar los primers válidos 
primers_validos (fw_primers, tga_primers, tag_primers, taa_primers)
#Guardar en un objeto
lista_primers <- primers_validos (fw_primers, tga_primers, tag_primers, taa_primers)
length (lista_primers)

print(lista_primers)

#Generar listas por de secuencias y nombres para el archivo FASTA

escribir_fasta <- function (lista_primers) {
  lista_secuencias <- list ()
  lista_nombres <- list ()   
  
  for (primer in lista_primers) {
    secuencias <- as.character (primer$secuencias) #agregué porque sino se imprime todo el DNAStringSet en el fasta T_T
    lista_secuencias <- c(lista_secuencias, secuencias)
    lista_nombres <- c(lista_nombres, primer$nombre)
  }
  
  if (length(lista_secuencias) > 0  & length(lista_nombres) > 0) {
    library(seqinr)
    write.fasta (sequences = lista_secuencias, 
                names = lista_nombres, 
                nbchar = 80, 
                file.out = "resultados/gallus_primer.fasta")
    print("Ver primers en carpeta de resultados")
  }
}

escribir_fasta (lista_primers)



