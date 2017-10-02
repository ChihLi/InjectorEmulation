############################################################
#              *** VALIDATION ***                          #
#    MADE BY CHIH-LI and SIMON ON 9/19/2017                #
#    NOTE: PLOT POD MODES OF CIRCUMFERENTIAL VELOCITY FLOW #
#          (FIGURE 6)                                      #
############################################################

r <- 3  ## circumferential velocity flow
load(paste0(output.folder.path, "/POD_expansion/", Response.ColumnIndex[r], "/mode.B.RData"))
load(paste0(output.folder.path, "/POD_expansion/", Response.ColumnIndex[r], "/mode.phi.RData"))

Coordinate <- rbind(a.coord, b.coord, d.coord, c.coord)
numCol <- 1000   #numCol colors in legend
colors <- matlab.like(numCol+1)

png(paste0(output.folder.path, "/POD_expansion/", Response.ColumnIndex[r], "/mode1.png"), width = 800, height = 400)
Range <- c(-6, -1)
y <- mode.phi[,1] * mode.B[1,1]
zcolor <- rep(0, nrow(mode.phi))
zcolor[y < Range[1]] <-  colors[1]
zcolor[y > Range[2]] <-  colors[numCol+1]
select.fg <- y > Range[1] & y < Range[2]
zcolor[select.fg] <- colors[(y[select.fg] - Range[1])/diff(Range)*numCol + 1] 
plot(Coordinate[,1], Coordinate[,2], col = zcolor, pch=15, cex=1.2,
     xlab = 'x (m)', ylab = 'y (m)', main = paste0("POD mode 1")) 
colorlegend(col = colors, zlim=range(Range))
dev.off()

png(paste0(output.folder.path, "/POD_expansion/", Response.ColumnIndex[r], "/mode2.png"), width = 800, height = 400)

Range <- c(-1, 1)
y <- mode.phi[,2] * mode.B[1,2]
zcolor <- rep(0, nrow(mode.phi))
zcolor[y < Range[1]] <-  colors[1]
zcolor[y > Range[2]] <-  colors[numCol+1]
select.fg <- y > Range[1] & y < Range[2]
zcolor[select.fg] <- colors[(y[select.fg] - Range[1])/diff(Range)*numCol + 1] 
plot(Coordinate[,1], Coordinate[,2], col = zcolor, pch=15, cex=1.2,
     xlab = 'x (m)', ylab = 'y (m)', main = paste0("POD mode 2")) 
colorlegend(col = colors, zlim=range(Range))
dev.off()