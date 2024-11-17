
library(Biostrings)

secun_prueba1 <- readDNAStringSet("extras/sequence2.fasta")
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


primer_fw(122, secun_prueba1)      
