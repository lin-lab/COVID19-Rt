library(ggplot2)

#define a helper function (borrowed from the "ez" package)
ezLev=function(x,new_order){
  for(i in rev(new_order)){
    x=relevel(x,ref=i)
  }
  return(x)
}

ggcorplot = function(data,var_text_size,cor_text_limits){
  # normalize data
  # for(i in 1:length(data)){
  #   data[,i]=(data[,i]-mean(data[,i]))/sd(data[,i])
  # }
  # obtain new data frame
  z=data.frame()
  i = 1
  j = i
  while(i<=length(data)){
    if(j>length(data)){
      i=i+1
      j=i
    }else{
      x = data[,i]
      y = data[,j]
      temp=as.data.frame(cbind(x,y))
      temp=cbind(temp,names(data)[i],names(data)[j])
      z=rbind(z,temp)
      j=j+1
    }
  }
  names(z)=c('x','y','x_lab','y_lab')
  z$x_lab = ezLev(factor(z$x_lab),names(data))
  z$y_lab = ezLev(factor(z$y_lab),names(data))
  z=z[z$x_lab!=z$y_lab,]
  #obtain correlation values
  z_cor = data.frame()
  i = 1
  j = i
  while(i<=length(data)){
    if(j>length(data)){
      i=i+1
      j=i
    }else{
      x = data[,i]
      y = data[,j]
      x_mid = min(x)+diff(range(x))/2
      y_mid = min(y)+diff(range(y))/2
      this_cor = cor(x,y)
      this_cor.test = cor.test(x,y)
      this_col = ifelse(this_cor.test$p.value<.05,'<.05','>.05')
      this_size = (this_cor)^2
      slope = lm(y~x)$coefficients[2]
      cor_text = paste(paste0("corr = ", format(round(this_cor, digits=2), nsmall = 2)),
                       paste0("slope = ", format(round(slope, digits=2), nsmall = 2)), sep = "\n")
      b=as.data.frame(cor_text)
      b=cbind(b,x_mid,y_mid,this_col,this_size,names(data)[j],names(data)[i])
      z_cor=rbind(z_cor,b)
      j=j+1
    }
  }
  names(z_cor)=c('cor','x_mid','y_mid','p','rsq','x_lab','y_lab')
  z_cor$x_lab = ezLev(factor(z_cor$x_lab),names(data))
  z_cor$y_lab = ezLev(factor(z_cor$y_lab),names(data))
  diag = z_cor[z_cor$x_lab==z_cor$y_lab,]
  z_cor=z_cor[z_cor$x_lab!=z_cor$y_lab,]
  f = facet_grid(y_lab~x_lab,scales='free')
  o = theme(
    panel.grid.minor = element_blank()
    ,panel.grid.major = element_blank()
    #,axis.ticks = element_blank()
    ,axis.text.y = element_text()
    ,axis.text.x = element_text()
    ,axis.title.y = element_text()
    ,axis.title.x = element_text()
    ,legend.position='none'
  )
  size_scale = scale_size(limits = c(0,1))
  return(
    ggplot(data = z, mapping = aes(x = x, y = y)) + 
      geom_point() + 
      geom_smooth(colour = 'blue', method = 'lm', fill = 'green') + 
      geom_text(geom_params = list(size=var_text_size)
                , data = diag
                , mapping = aes(
                  x=y_mid
                  , y=x_mid
                  , label=x_lab
                ))+geom_text(data = z_cor
                             , mapping = aes(
                               x=y_mid
                               , y=x_mid
                               , label=cor
                               , size = rsq
                               , colour = p
                             ))+
      scale_colour_manual(values=c('blue','blue'))+
      f+
      o+
      size_scale
  )
}

