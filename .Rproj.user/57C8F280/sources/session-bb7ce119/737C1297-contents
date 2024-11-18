detener <- 72
inicio <- 0
ultima <- 20
while (inicio <= detener & ultima < 92 ) {
  inicio <- inicio + 1
  
  primer_fd <- subseq(secun_prueba, start=inicio, end=ultima) 
  
  print(primer_fd)
  
  ultima <- ultima + 1
  
  ##### Eveluar condiciones de primer
  
}

primer_fd # Ultima secuencia que me dio


# INICIO



# FIN


### CONTENIDO DE CG
# 50 - 60 %
longt <- width(primer_fd)
longt

cont_cg <- letterFrequency(primer_fd, "CG")
cont_cg

porc_cg <- (cont_cg / longt) * 100
porc_cg

### No tener los PATRONES 
# GGG | CCC | GGT | ATT | GGA | TAA | TTA 

tripletes <- trinucleotideFrequency(primer_fd)

patrones_malos <- tripletes[c(22, 43, 44, 16, 41, 49, 61)]
patrones_malos

no_hay <- c(0,0,0,0,0,0,0)

all(patrones_malos == 0)

### TEMPERATURA
# Formula: TM = 4GC + 2AT
cont_cg <- letterFrequency(primer_fd, "CG")
cont_cg
cont_at <- letterFrequency(primer_fd, "AT")
cont_at

#  Tm entre los 55 - 65
temperatura <- (4*cont_cg) + (2*cont_at)
temperatura

cont_cg <- as.numeric(letterFrequency(primer_fd, "CG"))
cont_cg
cont_at <- as.numeric(letterFrequency(primer_fd, "AT"))
cont_at
