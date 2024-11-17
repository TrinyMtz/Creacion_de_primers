
library(Biostrings)

secun_prueba1 <- readDNAStringSet("Creacion_de_primers/extras/gallus.fasta")
secun_prueba1


pre_fw <- function(secun_prueba) {
  tipo_sec <- class(secun_prueba)
  tipo_sec[1]
  
  # Para evaluar tp de secuencia
  if (tipo_sec[1] == "DNAStringSet") {
    # Comprobar que no es degenerado
    no_degenerado <- alphabetFrequency(secun_prueba, baseOnly=TRUE)
    
    if (no_degenerado[5] == 0) {
      # Verificar el largo de la secuencia
      longt <- width(secun_prueba)
      if ( longt < 20000) {
      # Codon de inicio: TAC -> ATG
        codon_in <- vmatchPattern("ATG", secun_prueba)
        codon_in[[1]][1]
      } else { print("La capacidad maxima es de 20,000 nucleotidos")}}} else { print("Cambiar a DNA")}
  }
    
pre_fw(secun_prueba1)    
####################33
        
        
        
primer_fw <-function(inicio_codon, secun_prueba) {
  detener <- inicio_codon - 20
  inicio <- 0
  ultima <- 20
  while (inicio <= detener & ultima < inicio_codon ) {
          inicio <- inicio + 1  
          primer_fd <- subseq(secun_prueba, start=inicio, end=ultima) 
          ultima <- ultima + 1
          
          ##### Eveluar condiciones de primer
          # Patrones
          tripletes <- trinucleotideFrequency(primer_fd)
          patrones_malos <- tripletes[c(22, 43, 44, 16, 41, 49, 61)]
          no_hay <- c(0,0,0,0,0,0,0)
          comparacion <- all(patrones_malos == no_hay)
          
          if (comparacion == TRUE) { #Seguir evaluando el primer 
            
            # Poct de CG
            longt <- width(primer_fd)
            cont_cg <- letterFrequency(primer_fd, "CG")
            porc_cg <- (cont_cg / longt) * 100
            
            if (porc_cg < 60 & porc_cg > 49) {#Seguir evaluando
              #Temperatura: 55 - 65 °C
              cont_cg <- letterFrequency(primer_fd, "CG")
              cont_at <- letterFrequency(primer_fd, "AT")
              temperatura <- (4*cont_cg) + (2*cont_at)
              
              if (temperatura >54 & temperatura <66 ) {
                print(primer_fd)
                print(paste("Porcentaje de CG: ", porc_cg))
                print(paste("Tm: ", temperatura))
              }
              }
          }
          } 
  } 


primer_fw(32, secun_prueba1) ##poner el valor de  interger     

###############
###Codones de terminación: TGA,TAG,TAA

#####Primer con TGA
revertida<-reverse(secun_prueba1)
revertida

vmatchPattern("AGT",revertida)->tga
tga

primer_rev_ct1<-function(tga1,secrev) {
  detener <- tga1 - 20
  inicio <- 0
  ultima <- 20
  while (inicio <= detener & ultima < tga1 ) {
    inicio <- inicio + 1  
    primer_rv<- subseq(revertida, start=inicio, end=ultima) 
    ultima <- ultima + 1
    final<-complement(primer_rv)#la secuencia complentaria
    
    tripletes <- trinucleotideFrequency(final)#revisar patrones indeseados
    patrones_malos <- tripletes[c(22, 43, 44, 16, 41, 49, 61)]
    no_hay <- c(0,0,0,0,0,0,0)
    comparacion <- all(patrones_malos == no_hay)
    tolerancia<-sum(patrones_malos!=0)
    
    if (comparacion == TRUE | tolerancia == 1) { #Si no están esos patrones que siga evaluando
      longt <- width(primer_rv)#revisar %gc
      cont_cg <- letterFrequency(primer_rv, "CG")
      porc_cg <- (cont_cg / longt) * 100
      if (porc_cg < 60 & porc_cg > 49){ 
        cont_cg <- letterFrequency(primer_rv, "CG")
        cont_at <- letterFrequency(primer_rv, "AT")
        
        temperatura <- (4*cont_cg) + (2*cont_at) ##Tm 
        if (temperatura >54 & temperatura <66 ) {
          print(primer_rv)
          print(paste("Porcentaje de CG: ", porc_cg))
          print(paste("Tm: ", temperatura))
        }
      }
    }
  }
} 


primer_rev_ct1(114,revertida)#primers con TGA 

#####Primer con TAG

vmatchPattern("GAT",revertida)->tag
tag

primer_rev_ct2<-function(tag1,secrev) {
  detener <- tag1 - 17
  inicio <- 0
  ultima <- 17
  while (inicio <= detener & ultima < tag1 ) {
    inicio <- inicio + 1  
    primer_rv<- subseq(revertida, start=inicio, end=ultima) 
    ultima <- ultima + 1
    final<-complement(primer_rv)##secuencia complementaria
    
    tripletes <- trinucleotideFrequency(primer_rv)#revisar patrones indeseados
    patrones_malos <- tripletes[c(22, 43, 44, 16, 41, 49, 61)]
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
          print(primer_rv)
          print(paste("Porcentaje de CG: ", porc_cg))
          print(paste("Tm: ", temperatura))
        }
      }
    }
  } 
} 

primer_rev_ct2(53,revertida)#sujeto a que tenga este patrón



######Primers con TAA
vmatchPattern("AAT",revertida)->taa
taa

primer_rev_ct3<-function(taa1,secrev) {
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
          print(primer_rv)
          print(paste("Porcentaje de CG: ", porc_cg))
          print(paste("Tm: ", temperatura))
        }
      }
    }
  } 
}
primer_rev_ct3(37,revertida)


