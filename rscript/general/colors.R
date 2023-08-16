color = list()
color$parula = colorRampPalette( c('#352187','#0f5cdd','#1481d6','#06a4ca','#2eb7a4','#87bf77','#d1bb59','#fec832','#f9fb0e'))
color$jet = colorRampPalette( c('#201040','#302070','#104090','#0060b0','#00a0a0','#80c080','#d0c010','#f0c010','#d07010','#b02030','#701020'))
color$jetblue = colorRampPalette( c('#201040','#302070','#104090','#0060b0'))
color$jetred = colorRampPalette( c('#701020','#b02030','#d07010','#f0c010'))
color$red = colorRampPalette( c('#f7f7f7','#fddbc7','#f4a582','#d6604d','#b2182b','#67001f'))
color$blue = colorRampPalette( c('#f7f7f7','#d1e5f0','#92c5de','#4393c3','#2166ac','#053061'))
color$green = colorRampPalette( c('#f7f7f7','#80c080','#309050','#106030'))
color$orange = colorRampPalette( c('#f7f7f7','#f0c010','#d07010','#603010'))
color$bluered = colorRampPalette( c('#4393c3','#92c5de','#d1e5f0','#f7f7f7','#fddbc7','#f4a582','#d6604d'))
color$pinkmel = colorRampPalette( c('#FFAC81', '#FF928B', '#FEC3A6', '#EFE9AE', '#CDEAC0'))
color$pinkpa = colorRampPalette( c('#F08080', '#F4978E', '#F8AD9D', '#FBC4AB', '#FFDAB9'))
color$japa = colorRampPalette( c('#F7B267', '#F79D65', '#F4845F', '#F27059', '#F25C54'))
color$greenpa = colorRampPalette( c('#DAF3D7', '#E4FDE1', '#C6EDC3', '#A7DCA5', '#90CF8E'))
color$orpa = colorRampPalette( c('#FCAC5D', '#FCBC5D', '#FCCC5D', '#FCD45D', '#FCEC5D'))
color$adjust = function(cols, brightness = 1, alpha = 1) {
  r = rgb2hsv(col2rgb(cols))
  r[3,] = r[3,] * brightness
  if(alpha < 1) return(hsv(r[1,], r[2,], r[3,], alpha))
  else return(hsv(r[1,], r[2,], r[3,]))
}

color$printall = function() {
  n = length(color)
  plot()
  for(i in names(color)) {
    print(color[[i]](20))
  }
  dev.off()
}

