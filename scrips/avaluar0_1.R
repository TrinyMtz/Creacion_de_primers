# PRIMER DE 20 b

detener <- 72
inicio <- 0
ultima <- 20

while (inicio <= detener & ultima < 92 ) {
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

