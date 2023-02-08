# Draw Venn-diagram like figure
# Should be in directory with folders 'Analytics', 'Description' etc. Usually called 'full_model'


library(magick)

# Load images if not loaded already
inputs=c("AE2","PIS","SMR00","SMR01","SMR01E","SMR04")
input_names=c("A&E","Prescriptions","Outpatients",
              "Hosp. inp./\nday case","Other",
              "MH inp./\nday case")

# Order, starting from rightmost and going clockwise
p_order=c(6,4,5,2,1,3)
inputs=inputs[p_order]
input_names=input_names[p_order]

if (!exists("panel6a")) {
  
  # General images
  panel0a = image_read_pdf("Description/age_all.pdf")
  panel0b = image_read_pdf("Description/simd_all.pdf")
  
  for (i in 1:6) {
    l1=list.files("Description",pattern=paste0("age_",inputs[i],"_."),full=T); f1=grep("pdf",l1,val=T)
    l2=list.files("Description",pattern=paste0("simd_",inputs[i],"_."),full=T); f2=grep("pdf",l2,val=T)
    panela=image_read_pdf(f1)
    panelb=image_read_pdf(f2)
    assign(paste0("panel",i,"a"),panela)
    assign(paste0("panel",i,"b"),panelb)
    print(paste0("Loaded images for ",input_names[i]))
  }    
}

# Read numbers
fnum=readLines("Description/misc.txt")
Rs=rep(0,length(inputs))
Is=rep(0,length(inputs))
for (i in 1:length(inputs)) {
  w=which(fnum==inputs[i])
  Rs[i]=fnum[w+1]
  Is[i]=fnum[w+2]
}
w1=which(fnum=="Total recorded episodes of any type")
Rtot=fnum[w1+1]
w2=which(fnum=="Total patients")
Itot=fnum[w2+1]

# Systemwatch- special case
w3=which(fnum=="SystemWatch")
Rsw=fnum[w3+1]
Isw=fnum[w3+2]


for (i in 1:length(inputs)) {
  Rs[i]=format(as.numeric(Rs[i]),big.mark=",")
  Is[i]=format(as.numeric(Is[i]),big.mark=",")
}
Rtot=format(as.numeric(Rtot),big.mark=","); Itot=format(as.numeric(Itot),big.mark=",")
Rsw=format(as.numeric(Rsw),big.mark=","); Isw=format(as.numeric(Isw),big.mark=",")



# Parameters
hcent = 0.2 # Horizontal offset for central plots
hsc = 0.2 # Scale for central plots

ocent=0.8 # offset from centre for radial plots
osc = 0.13 # Scale for radial plots

pdf("../../Diagrams/venn.pdf",width=4,height=4)
par(mar=rep(0.1,4))
plot(0,type="n",xlim=c(-1,1),ylim=c(-1,1),xaxt="n",yaxt="n",ann=F,bty="n")

# Circle
tt=seq(0,2*pi,length=200); sc=0.5
lines(sc*sin(tt),sc*cos(tt),lwd=2)

# Central plots
rasterImage(panel0a,-hcent-hsc,-hsc,-hcent+hsc,hsc)
rasterImage(panel0b,hcent-hsc,-hsc,hcent+hsc,hsc)

# Text
text(0,hsc + (sc-hsc)/2,"Total",adj=c(0.5,0.5))
text(0,-hsc - (sc-hsc)/2,paste0("Records (R) = ",Rtot," \n Individuals (I) = ",Itot),
     cex=0.5,adj=c(0.5,0.5))

# Other plots
pnames=""

for (i in 1:6) {

  angle=2*(i+0.5)*pi/6
  xcentre=ocent*sin(angle)
  ycentre=ocent*cos(angle)
  
  # Get figures
  panela=get(paste0("panel",i,"a"))
  panelb=get(paste0("panel",i,"b"))
  
  # Plots
  rasterImage(panela,xcentre-2*osc,ycentre-osc,xcentre,ycentre+osc)
  rasterImage(panelb,xcentre,ycentre-osc,xcentre+2*osc,ycentre+osc)
  
  # Text
  if (input_names[i]!="Other") xtext=paste0("R = ",Rs[i]," \n I = ",Is[i]) else {
    xtext=paste0("R = ",Rs[i]," (Ger.); ",Rsw," (Sys.)\nI = ",Is[i]," (Ger.); ",Isw," (Sys.)")
  }
  text(xcentre,ycentre + osc*1.5,input_names[i],adj=c(0.5,0.5),cex=0.5)
  text(xcentre,ycentre - osc*1.5,xtext,cex=0.5,adj=c(0.5,0.5))

}
dev.off()