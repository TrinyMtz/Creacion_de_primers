library(Biostrings)

# Lectura de la secuencia
secun_prueba1 <- readDNAStringSet("Creacion_de_primers/extras/musmusculusBAX.fasta")
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
        codon_in <- vmatchPattern("ATG", secun_prueba)
        codon_in[[1]][1]
        print(codon_in)
      } else { print("La capacidad maxima es de 20,000 nucleotidos")}}} else { print("Cambiar a DNA")}
}

pre_fw(secun_prueba1) # Prueba

####################

#  GENERADOR DE PRIMERS FORWARD

primer_fw <-function(inicio_codon, secun_prueba, ultima) {
  detener <- inicio_codon - ultima
  inicio <- 0
  while (inicio <= detener & ultima < inicio_codon ) {
    inicio <- inicio + 1  
    primer_fd <- subseq(secun_prueba, start=inicio, end=ultima) 
    ultima <- ultima + 1
    
    final<-complement(primer_fd)
    ##### Eveluar condiciones de primer
    # Patrones
    tripletes <- trinucleotideFrequency(final)
    patrones_malos <- tripletes[c(22, 43, 44, 16, 41, 49, 61)]
    no_hay <- c(0,0,0,0,0,0,0)
    comparacion <- all(patrones_malos == no_hay)
    tolerancia<-sum(patrones_malos!=0)
    if (comparacion == TRUE |tolerancia==1) { #Seguir evaluando el primer 
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
          
          print("Forward")
          print(final)
          print(paste("Porcentaje de CG: ", porc_cg))
          print(paste("Tm: ", temperatura))
          
        }
      }
    }
  } 
} 


primer_fw(63, secun_prueba1, 21) # Colocar el valor del start , secuencia, long del primer. 
#Hasta 150 encuentra 2, pero se pasa por más de la mitad del límite

############################
####################     PRIMERS REVERSE      ############################
############################

### Codones de terminación: TGA,TAG,TAA

####################
# Primer con TGA
revertida <-reverse(secun_prueba1)
revertida

vmatchPattern("AGT",revertida)->tga
tga

primer_rev_ct1<-function(tga1,secrev,ultima) {
  detener <- tga1 - ultima
  inicio <- 0
  while (inicio <= detener & ultima < tga1 ) {
    inicio <- inicio + 1  
    primer_rv<- subseq(revertida, start=inicio, end=ultima) 
    ultima <- ultima + 1
    final<-complement(primer_rv)#la secuencia complentaria
    
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
          print("reverse tga")
          print(primer_rv)
          print(paste("Porcentaje de CG: ", porc_cg))
          print(paste("Tm: ", temperatura))
        }
      }
    }
  }
} 

primer_rev_ct1(100,revertida, 20) # Primers con TGA, puedes, jugar con los valores de 
# Start del objeto TGA, considerando tu región codificante (CDS)

####################
# Primer con TAG

vmatchPattern("GAT",revertida)->tag
tag

primer_rev_ct2<-function(tag1,secrev,ultima) {
  detener <- tag1 - ultima
  inicio <- 0
  while (inicio <= detener & ultima < tag1 ) {
    inicio <- inicio + 1  
    primer_rv<- subseq(revertida, start=inicio, end=ultima) 
    ultima <- ultima + 1
    final<-complement(primer_rv)##secuencia complementaria
    
    tripletes <- trinucleotideFrequency(primer_rv)#revisar patrones indeseados
    patrones_malos <- tripletes[c(22, 43, 44, 16, 41, 49, 61)]
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
          print("reverse tag")
          print(primer_rv)
          print(paste("Porcentaje de CG: ", porc_cg))
          print(paste("Tm: ", temperatura))
        }
      }
    }
  } 
} 

primer_rev_ct2(387,revertida, 20)#reverse TAG,sujeto a que tenga este patrón,puedes
#jugar con los valores de start del objeto tga, considerando tu región codificante (CDS)


####################
# Primers con TAA
vmatchPattern("AAT",revertida)->taa
taa

primer_rev_ct3<-function(taa1,secrev,ultima) {
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
          print("reverse taa")
          print(primer_rv)
          print(paste("Porcentaje de CG: ", porc_cg))
          print(paste("Tm: ", temperatura))
        }
      }
    }
  } 
}

primer_rev_ct3(18,revertida, 18)#reverse TAG,sujeto a que tenga este patrón,puedes
#jugar con los valores de start del objeto tga, considerando tu región codificante (CDS)


#####Con esta secuencia solo sale un reverse, los forward no



#####Prueba 2


# Lectura de la secuencia
secun_prueba1 <- readDNAStringSet("Creacion_de_primers/extras/frogESR.fasta")
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
        codon_in <- vmatchPattern("ATG", secun_prueba)
        codon_in[[1]][1]
      } else { print("La capacidad maxima es de 20,000 nucleotidos")}}} else { print("Cambiar a DNA")}
}

pre_fw(secun_prueba1) # Prueba

####################

#  GENERADOR DE PRIMERS FORWARD

primer_fw <-function(inicio_codon, secun_prueba, ultima) {
  detener <- inicio_codon - ultima
  inicio <- 0
  while (inicio <= detener & ultima < inicio_codon ) {
    inicio <- inicio + 1  
    primer_fd <- subseq(secun_prueba, start=inicio, end=ultima) 
    ultima <- ultima + 1
    
    final<-complement(primer_fd)
    ##### Eveluar condiciones de primer
    # Patrones
    tripletes <- trinucleotideFrequency(final)
    patrones_malos <- tripletes[c(22, 43, 44, 16, 41, 49, 61)]
    no_hay <- c(0,0,0,0,0,0,0)
    comparacion <- all(patrones_malos == no_hay)
    tolerancia<-sum(patrones_malos!=0)
    if (comparacion == TRUE |tolerancia==1) { #Seguir evaluando el primer 
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
          
          print("Forward")
          print(final)
          print(paste("Porcentaje de CG: ", porc_cg))
          print(paste("Tm: ", temperatura))
          
        }
      }
    }
  } 
} 


primer_fw(137, secun_prueba1, 20) # Colocar el valor del start , secuencia, long del primer. 
#Hasta 150 encuentra 2, pero se pasa por más de la mitad del límite

############################
####################     PRIMERS REVERSE      ############################
############################

### Codones de terminación: TGA,TAG,TAA

####################
# Primer con TGA
revertida <-reverse(secun_prueba1)
revertida

vmatchPattern("AGT",revertida)->tga
tga

primer_rev_ct1<-function(tga1,secrev,ultima) {
  detener <- tga1 - ultima
  inicio <- 0
  while (inicio <= detener & ultima < tga1 ) {
    inicio <- inicio + 1  
    primer_rv<- subseq(revertida, start=inicio, end=ultima) 
    ultima <- ultima + 1
    final<-complement(primer_rv)#la secuencia complentaria
    
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
          print("reverse tga")
          print(primer_rv)
          print(paste("Porcentaje de CG: ", porc_cg))
          print(paste("Tm: ", temperatura))
        }
      }
    }
  }
} 

primer_rev_ct1(29,revertida, 20) # Primers con TGA, puedes, jugar con los valores de 
# Start del objeto TGA, considerando tu región codificante (CDS)

####################
# Primer con TAG

vmatchPattern("GAT",revertida)->tag
tag

primer_rev_ct2<-function(tag1,secrev,ultima) {
  detener <- tag1 - ultima
  inicio <- 0
  while (inicio <= detener & ultima < tag1 ) {
    inicio <- inicio + 1  
    primer_rv<- subseq(revertida, start=inicio, end=ultima) 
    ultima <- ultima + 1
    final<-complement(primer_rv)##secuencia complementaria
    
    tripletes <- trinucleotideFrequency(primer_rv)#revisar patrones indeseados
    patrones_malos <- tripletes[c(22, 43, 44, 16, 41, 49, 61)]
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
          print("reverse tag")
          print(primer_rv)
          print(paste("Porcentaje de CG: ", porc_cg))
          print(paste("Tm: ", temperatura))
        }
      }
    }
  } 
} 

primer_rev_ct2(106,revertida, 20)#reverse TAG,sujeto a que tenga este patrón,puedes
#jugar con los valores de start del objeto tga, considerando tu región codificante (CDS)


####################
# Primers con TAA
vmatchPattern("AAT",revertida)->taa
taa

primer_rev_ct3<-function(taa1,secrev,ultima) {
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
          print("reverse taa")
          print(primer_rv)
          print(paste("Porcentaje de CG: ", porc_cg))
          print(paste("Tm: ", temperatura))
        }
      }
    }
  } 
}

primer_rev_ct3(100,revertida, 24)#reverse TAG,sujeto a que tenga este patrón,puedes
#jugar con los valores de start del objeto tga, considerando tu región codificante (CDS)

##En esta secuencia por ejemplo, en el extremo 5´-3´solo hay 9 bases, tiene que pasarse,
#es desición del ususario cuanto

#prueba 3

# Lectura de la secuencia
secun_prueba1 <- readDNAStringSet("Creacion_de_primers/extras/cabraaopae.fasta")
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
        codon_in <- vmatchPattern("ATG", secun_prueba)
        codon_in[[1]][1]
        print(codon_in)
      } else { print("La capacidad maxima es de 20,000 nucleotidos")}}} else { print("Cambiar a DNA")}
}

pre_fw(secun_prueba1) #Identifica un atg bastante lejos del inicio de la región codificante

####################

#  GENERADOR DE PRIMERS FORWARD

primer_fw <-function(inicio_codon, secun_prueba, ultima) {
  detener <- inicio_codon - ultima
  inicio <- 0
  while (inicio <= detener & ultima < inicio_codon ) {
    inicio <- inicio + 1  
    primer_fd <- subseq(secun_prueba, start=inicio, end=ultima) 
    ultima <- ultima + 1
    
    final<-complement(primer_fd)
    ##### Eveluar condiciones de primer
    # Patrones
    tripletes <- trinucleotideFrequency(final)
    patrones_malos <- tripletes[c(22, 43, 44, 16, 41, 49, 61)]
    no_hay <- c(0,0,0,0,0,0,0)
    comparacion <- all(patrones_malos == no_hay)
    tolerancia<-sum(patrones_malos!=0)
    if (comparacion == TRUE |tolerancia==1) { #Seguir evaluando el primer 
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
          
          print("Forward")
          print(final)
          print(paste("Porcentaje de CG: ", porc_cg))
          print(paste("Tm: ", temperatura))
          
        }
      }
    }
  } 
} 


primer_fw(100, secun_prueba1,19) # Colocar el valor del start , secuencia, long del primer. 
#Hasta 150 encuentra 2, pero se pasa por más de la mitad del límite

############################
####################     PRIMERS REVERSE      ############################
############################

### Codones de terminación: TGA,TAG,TAA

####################
# Primer con TGA
revertida <-reverse(secun_prueba1)
revertida

vmatchPattern("AGT",revertida)->tga
tga

primer_rev_ct1<-function(tga1,secrev,ultima) {
  detener <- tga1 - ultima
  inicio <- 0
  while (inicio <= detener & ultima < tga1 ) {
    inicio <- inicio + 1  
    primer_rv<- subseq(revertida, start=inicio, end=ultima) 
    ultima <- ultima + 1
    final<-complement(primer_rv)#la secuencia complentaria
    
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
          print("reverse tga")
          print(primer_rv)
          print(paste("Porcentaje de CG: ", porc_cg))
          print(paste("Tm: ", temperatura))
        }
      }
    }
  }
} 

primer_rev_ct1(142,revertida, 20) # Primers con TGA, puedes, jugar con los valores de 
# Start del objeto TGA, considerando tu región codificante (CDS)

####################
# Primer con TAG

vmatchPattern("GAT",revertida)->tag
tag

primer_rev_ct2<-function(tag1,secrev,ultima) {
  detener <- tag1 - ultima
  inicio <- 0
  while (inicio <= detener & ultima < tag1 ) {
    inicio <- inicio + 1  
    primer_rv<- subseq(revertida, start=inicio, end=ultima) 
    ultima <- ultima + 1
    final<-complement(primer_rv)##secuencia complementaria
    
    tripletes <- trinucleotideFrequency(primer_rv)#revisar patrones indeseados
    patrones_malos <- tripletes[c(22, 43, 44, 16, 41, 49, 61)]
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
          print("reverse tag")
          print(primer_rv)
          print(paste("Porcentaje de CG: ", porc_cg))
          print(paste("Tm: ", temperatura))
        }
      }
    }
  } 
} 

primer_rev_ct2(31,revertida, 20)#reverse TAG,sujeto a que tenga este patrón,puedes
#jugar con los valores de start del objeto tga, considerando tu región codificante (CDS)


####################
# Primers con TAA
vmatchPattern("AAT",revertida)->taa
taa

primer_rev_ct3<-function(taa1,secrev,ultima) {
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
          print("reverse taa")
          print(primer_rv)
          print(paste("Porcentaje de CG: ", porc_cg))
          print(paste("Tm: ", temperatura))
        }
      }
    }
  } 
}

primer_rev_ct3(96,revertida,17)#reverse TAG,sujeto a que tenga este patrón,puedes
#jugar con los valores de start del objeto tga, considerando tu región codificante (CDS)



#########Prueba 4
# Lectura de la secuencia
secun_prueba1 <- readDNAStringSet("Creacion_de_primers/extras/aktsusscrofa.fasta")
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
        codon_in <- vmatchPattern("ATG", secun_prueba)
        codon_in[[1]][1]
        print(codon_in)
      } else { print("La capacidad maxima es de 20,000 nucleotidos")}}} else { print("Cambiar a DNA")}
}

pre_fw(secun_prueba1) #Si el límite por el codon de inicio es muy corto, explorar otras posiciones

####################

#  GENERADOR DE PRIMERS FORWARD

primer_fw <-function(inicio_codon, secun_prueba, ultima) {
  detener <- inicio_codon - ultima
  inicio <- 0
  while (inicio <= detener & ultima < inicio_codon ) {
    inicio <- inicio + 1  
    primer_fd <- subseq(secun_prueba, start=inicio, end=ultima) 
    ultima <- ultima + 1
    
    final<-complement(primer_fd)
    ##### Eveluar condiciones de primer
    # Patrones
    tripletes <- trinucleotideFrequency(final)
    patrones_malos <- tripletes[c(22, 43, 44, 16, 41, 49, 61)]
    no_hay <- c(0,0,0,0,0,0,0)
    comparacion <- all(patrones_malos == no_hay)
    tolerancia<-sum(patrones_malos!=0)
    if (comparacion == TRUE |tolerancia==1) { #Seguir evaluando el primer 
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
          
          print("Forward")
          print(final)
          print(paste("Porcentaje de CG: ", porc_cg))
          print(paste("Tm: ", temperatura))
          
        }
      }
    }
  } 
} 


primer_fw(55, secun_prueba1, 20) # Colocar el valor del start , secuencia, long del primer. 
#Hasta 150 encuentra 2, pero se pasa por más de la mitad del límite

#Explorar otras posiciones del atg, da resultados

############################
####################     PRIMERS REVERSE      ############################
############################

### Codones de terminación: TGA,TAG,TAA

####################
# Primer con TGA
revertida <-reverse(secun_prueba1)
revertida

vmatchPattern("AGT",revertida)->tga
tga

primer_rev_ct1<-function(tga1,secrev,ultima) {
  detener <- tga1 - ultima
  inicio <- 0
  while (inicio <= detener & ultima < tga1 ) {
    inicio <- inicio + 1  
    primer_rv<- subseq(revertida, start=inicio, end=ultima) 
    ultima <- ultima + 1
    final<-complement(primer_rv)#la secuencia complentaria
    
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
          print("reverse tga")
          print(primer_rv)
          print(paste("Porcentaje de CG: ", porc_cg))
          print(paste("Tm: ", temperatura))
        }
      }
    }
  }
} 

primer_rev_ct1(383,revertida, 20) # Primers con TGA, puedes, jugar con los valores de 
# Start del objeto TGA, considerando tu región codificante (CDS)

####################
# Primer con TAG

vmatchPattern("GAT",revertida)->tag
tag

primer_rev_ct2<-function(tag1,secrev,ultima) {
  detener <- tag1 - ultima
  inicio <- 0
  while (inicio <= detener & ultima < tag1 ) {
    inicio <- inicio + 1  
    primer_rv<- subseq(revertida, start=inicio, end=ultima) 
    ultima <- ultima + 1
    final<-complement(primer_rv)##secuencia complementaria
    
    tripletes <- trinucleotideFrequency(primer_rv)#revisar patrones indeseados
    patrones_malos <- tripletes[c(22, 43, 44, 16, 41, 49, 61)]
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
          print("reverse tag")
          print(primer_rv)
          print(paste("Porcentaje de CG: ", porc_cg))
          print(paste("Tm: ", temperatura))
        }
      }
    }
  } 
} 

primer_rev_ct2(509,revertida, 20)#reverse TAG,sujeto a que tenga este patrón,puedes
#jugar con los valores de start del objeto tga, considerando tu región codificante (CDS)


####################
# Primers con TAA
vmatchPattern("AAT",revertida)->taa
taa

primer_rev_ct3<-function(taa1,secrev,ultima) {
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
          print("reverse taa")
          print(primer_rv)
          print(paste("Porcentaje de CG: ", porc_cg))
          print(paste("Tm: ", temperatura))
        }
      }
    }
  } 
}

primer_rev_ct3(100,revertida, 24)#reverse TAG,sujeto a que tenga este patrón,puedes
#jugar con los valores de start del objeto tga, considerando tu región codificante (CDS)