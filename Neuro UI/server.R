

library(shiny)
library(astsa)
library(igraph)
library(eegkit)
library(magick)
library(magrittr)
library(animation)
#Return the spectral matrices of the band chosen by the user
spec_matrix=function(data,bande){
  
  # Getting the frequencies from mvspec 
  u=mvspec(data,plot = F)
  X=u$freq
  
  #Getting all the spectral matrices for different frequencies 
  M=u$fxx
  
  #Delta waves
  if (bande==1){
    
    spec=M[,,X>=.1/128 & X<=3/128]
    return(spec)
  }
  #Theta waves
  else if (bande==2){
    
    spec=M[,,X>=4/128 & X<=7/128]
    return(spec)
  }
  #Alpha waves
  else if (bande==3){
    
    spec=M[,,X>=8/128 & X<=15/128]
    return(spec)
  }
  #Beta waves
  else if (bande==4){
    spec=M[,,X>=16/128 & X<=30/128]
    return(spec)
  }
  #Gamma Waves
  else {
    spec=M[,,X>30/128 & X<=64/128]
    return(spec)
  }
  
}

#Computing the coherence matrix from the spectral matrix 
coh_matrix=function(spect){
    n=length(spect[,1])
    cohh=matrix(rep(0),nrow = n,ncol = n)
    d=Re(diag(spect))
    dd=tcrossprod(d)
    cohh=Mod(spect)**2/dd
    return(cohh)
}

#Compute the mean of the coherence matrices
moy_coh=function(all_spect_matrix){
    return(coh_matrix(apply(all_spect_matrix,c(1,2),mean)))
}


#Computing the adjacency matrix from the coherence matrix and the threshold chosen by the user 
adj_matrix=function(thresh,M){
    n=length(M[,1])
    adjm=matrix(as.numeric(M>thresh),nrow = n,ncol = n)
    return(adjm)
}

connectivity_power=function(graph,coherence_matrix){
  
  #getting the vertices of each edge in the graph
  d=ends(graph,E(graph))
  
  #Getting the corresponding coherence values
  cv=coherence_matrix[d]
  n=length(cv)
  #Setting the power of connectivity for each edge
  pow=rep(1,n)
  pow[cv>.6 & cv<=.7]=2
  pow[cv>.7 & cv<=.8]=4
  pow[cv>.8 & cv<=.9]=6
  pow[cv>.9 & cv<=1]=8
  g <- graph %>%
    set_edge_attr("weight", value = pow) %>%
    set_edge_attr("color", value = "red")
  g
  return(g)
}

#Creating animation from the interval given by the user.
anim=function(eeg_data_c,eeg_data_o,interval,thresh,bande,name){
  #Getting the name of the electrodes of interest.
  myelectrodes <- colnames(eeg_data_c)
  data("eegcoord")
  
  #Getting the coordinates of the electrodes of interest.
  two_Dcoord=as.matrix(eegcoord[rownames(eegcoord)%in%myelectrodes,c(4,5)])
  
  #Sorting the data to get the correct placement of electrodes
  elec_names=rownames(two_Dcoord)
  eeg_data_c=eeg_data_c[,elec_names]
  eeg_data_o=eeg_data_o[,elec_names]
  
  n=min(length(eeg_data_c[,1]),length(eeg_data_o[,1]))
  k=interval*128
  d=seq(1,n,k)
  
  for (i in 1:(length(d)-1)){
    #Spectral matrix for eyes closed
    specu=spec_matrix(eeg_data_c[d[i]:d[i+1],],bande)
    #Coherence Matrix
    cohu=moy_coh(specu)
    #Adjacency matrix for eyes closed
    adjmu=adj_matrix(thresh,cohu)
    #Graph from the adjacency matrix for eyes closed
    g1 <- graph_from_adjacency_matrix(adjmu,mode = "undirected",diag=F)
    g1 <-connectivity_power(g1,cohu)
    #Spectral matrix for eyes open
    specv=spec_matrix(eeg_data_o[d[i]:d[i+1],],bande)
    #Coherence Matrix
    cohv=moy_coh(specv)
    #Adjacency matrix for eyes open
    adjmv=adj_matrix(thresh,cohv)
    #Graph from the adjacency matrix for eyes open
    g2 <- graph_from_adjacency_matrix(adjmv,mode = "undirected",diag=F)
    g2 <-connectivity_power(g2,cohv)
    par(mfrow=c(1,2))
    plot(g1,layout=two_Dcoord,vertex.color="pink",
         edge.color="black",edge.width=E(g1)$weight,
         vertex.label=elec_names,main="Closed")
    plot(g2,layout=two_Dcoord,vertex.color="pink",edge.width=E(g2)$weight,
         edge.color="red",vertex.label=elec_names,main="Open")
  }
  
}



# Define server to get the coherence and adjacency matrices and plot the graphs
shinyServer(function(input, output) {
  #Getting the data for eyes closed 
  mind=reactive({ 
    
    f1=input$EEG_Closed
    
    #Loading the Data 
    EEG_closed=read.table(f1$datapath,sep=input$data_type,header =T)
    #Getting rid of the first column of the dataframe (Time).
    EEG_closed=EEG_closed[-1]
    
    #Getting the data for eyes open
    f2=input$EEG_Open
    
    #Loading the Data 
    EEG_open=read.table(f2$datapath,sep=input$data_type,header = T)
    EEG_open=EEG_open[-1]
    dir=getwd()
    if (grepl("/www",dir,fixed=T )){setwd(dir)}
    else{setwd(paste(dir,"/www",sep=""))}
    
    saveVideo(
      anim(EEG_closed,EEG_open,input$interval,input$thresh,input$band,"open")
      ,
      video.name = paste(input$name,".mp4",sep=""), 
      img.name = "J1",
      interval = 1, 
      ani.width = 1080, 
      ani.height = 640,
      outdir =getwd()
    )
  })
  
  
  
  #Plotting both graphs eyes closed on the left   and eyes open on the right.
  output$distPlot <- renderUI({
    validate(
      need(input$EEG_Closed != "", "Please select The Eye Open Data:")
    )
    validate(
      need(input$EEG_Open != "", "Please select the Eye Closed Data:")
    )
    
    if (input$go>0) { 
      isolate(mind())
      isolate(h6("Visualization of Spikes", br(),
         tags$video(src=paste(input$name,".mp4",sep=""),type="video/mp4",
                    width="720px",heigth="640px",controls="controls")))
    } 
    
    
  })
  
})
