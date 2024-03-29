// Generated by CoffeeScript 1.12.6
var iheatmap;

iheatmap = function(widgetdiv, data, chartOpts) {
  var axispos, cells, chartdivid, colors, flip_vert_slice, formatX, formatY, g_heatmap, g_horslice, g_verslice, hbot, height, horcurve, horslice, htop, linecolor, linewidth, margin, myheatmap, nullcolor, nxticks, nyticks, nzticks, plotHor, plotVer, rectcolor, ref, ref1, ref10, ref11, ref12, ref13, ref14, ref15, ref16, ref17, ref18, ref19, ref2, ref20, ref21, ref22, ref23, ref24, ref25, ref26, ref27, ref28, ref29, ref3, ref30, ref31, ref4, ref5, ref6, ref7, ref8, ref9, svg, title, titlepos, ver_opts, vercurve, verslice, widgetdivid, width, wleft, wright, xdif, xlab, xlim, xticks, ydif, ylab, ylim, yticks, z_transpose, zlab, zlim, zthresh, zticks;
  height = (ref = chartOpts != null ? chartOpts.height : void 0) != null ? ref : 800;
  width = (ref1 = chartOpts != null ? chartOpts.width : void 0) != null ? ref1 : 800;
  htop = (ref2 = chartOpts != null ? chartOpts.htop : void 0) != null ? ref2 : height / 2;
  wleft = (ref3 = chartOpts != null ? chartOpts.wleft : void 0) != null ? ref3 : width / 2;
  margin = (ref4 = chartOpts != null ? chartOpts.margin : void 0) != null ? ref4 : {
    left: 60,
    top: 40,
    right: 40,
    bottom: 40,
    inner: 0
  };
  axispos = (ref5 = chartOpts != null ? chartOpts.axispos : void 0) != null ? ref5 : {
    xtitle: 25,
    ytitle: 30,
    xlabel: 5,
    ylabel: 5
  };
  titlepos = (ref6 = chartOpts != null ? chartOpts.titlepos : void 0) != null ? ref6 : 20;
  rectcolor = (ref7 = chartOpts != null ? chartOpts.rectcolor : void 0) != null ? ref7 : "#E6E6E6";
  nullcolor = (ref8 = chartOpts != null ? chartOpts.nullcolor : void 0) != null ? ref8 : "#E6E6E6";
  linecolor = (ref9 = chartOpts != null ? chartOpts.linecolor : void 0) != null ? ref9 : "slateblue";
  linewidth = (ref10 = chartOpts != null ? chartOpts.linewidth : void 0) != null ? ref10 : 2;
  xlim = (ref11 = chartOpts != null ? chartOpts.xlim : void 0) != null ? ref11 : null;
  ylim = (ref12 = chartOpts != null ? chartOpts.ylim : void 0) != null ? ref12 : null;
  nxticks = (ref13 = chartOpts != null ? chartOpts.nxticks : void 0) != null ? ref13 : 5;
  xticks = (ref14 = chartOpts != null ? chartOpts.xticks : void 0) != null ? ref14 : null;
  nyticks = (ref15 = chartOpts != null ? chartOpts.nyticks : void 0) != null ? ref15 : 5;
  yticks = (ref16 = chartOpts != null ? chartOpts.yticks : void 0) != null ? ref16 : null;
  nzticks = (ref17 = chartOpts != null ? chartOpts.nzticks : void 0) != null ? ref17 : 5;
  zticks = (ref18 = chartOpts != null ? chartOpts.zticks : void 0) != null ? ref18 : null;
  title = (ref19 = chartOpts != null ? chartOpts.title : void 0) != null ? ref19 : "";
  xlab = (ref20 = chartOpts != null ? chartOpts.xlab : void 0) != null ? ref20 : "X";
  ylab = (ref21 = chartOpts != null ? chartOpts.ylab : void 0) != null ? ref21 : "Y";
  zlab = (ref22 = chartOpts != null ? chartOpts.zlab : void 0) != null ? ref22 : "Z";
  zthresh = (ref23 = chartOpts != null ? chartOpts.zthresh : void 0) != null ? ref23 : null;
  zlim = (ref24 = chartOpts != null ? chartOpts.zlim : void 0) != null ? ref24 : [-d3panels.matrixMaxAbs(data.z), 0, d3panels.matrixMaxAbs(data.z)];
  colors = (ref25 = chartOpts != null ? chartOpts.colors : void 0) != null ? ref25 : ["slateblue", "white", "crimson"];
  flip_vert_slice = (ref26 = chartOpts != null ? chartOpts.flip_vert_slice : void 0) != null ? ref26 : false;
  chartdivid = (ref27 = chartOpts != null ? chartOpts.chartdivid : void 0) != null ? ref27 : 'chart';
  widgetdivid = d3.select(widgetdiv).attr('id');
  margin = d3panels.check_listarg_v_default(margin, {
    left: 60,
    top: 40,
    right: 40,
    bottom: 40,
    inner: 0
  });
  axispos = d3panels.check_listarg_v_default(axispos, {
    xtitle: 25,
    ytitle: 30,
    xlabel: 5,
    ylabel: 5
  });
  hbot = height - htop;
  wright = width - wleft;
  svg = d3.select(widgetdiv).select("svg");
  if (xlim == null) {
    xlim = d3.extent(data.x);
    xdif = (data.x[1] - data.x[0]) / 2;
    xlim[0] -= xdif;
    xlim[1] += xdif;
  }
  if (ylim == null) {
    ylim = d3.extent(data.y);
    ydif = (data.y[1] - data.y[0]) / 2;
    ylim[0] -= ydif;
    ylim[1] += ydif;
  }
  z_transpose = d3panels.transpose(data.z);
  myheatmap = d3panels.heatmap({
    width: wleft,
    height: htop,
    margin: margin,
    axispos: axispos,
    titlepos: titlepos,
    rectcolor: rectcolor,
    xlim: xlim,
    ylim: ylim,
    nxticks: nxticks,
    xticks: xticks,
    nyticks: nyticks,
    yticks: yticks,
    xlab: xlab,
    ylab: ylab,
    zlim: zlim,
    zthresh: zthresh,
    colors: colors,
    nullcolor: nullcolor,
    tipclass: widgetdivid
  });
  horslice = d3panels.panelframe({
    width: wleft,
    height: hbot,
    margin: margin,
    axispos: axispos,
    titlepos: titlepos,
    rectcolor: rectcolor,
    xlim: xlim,
    ylim: d3.extent(zlim),
    nxticks: nxticks,
    xticks: xticks,
    nyticks: nzticks,
    yticks: zticks,
    xlab: xlab,
    ylab: zlab
  });
  ver_opts = {
    width: wright,
    height: htop,
    margin: margin,
    axispos: axispos,
    titlepos: titlepos,
    rectcolor: rectcolor,
    xlim: ylim,
    ylim: d3.extent(zlim),
    nxticks: nyticks,
    xticks: yticks,
    nyticks: nzticks,
    yticks: zticks,
    xlab: ylab,
    ylab: zlab
  };
  if (flip_vert_slice) {
    ref28 = [ver_opts.ylab, ver_opts.xlab], ver_opts.xlab = ref28[0], ver_opts.ylab = ref28[1];
    ref29 = [ver_opts.ylim, ver_opts.xlim], ver_opts.xlim = ref29[0], ver_opts.ylim = ref29[1];
    ref30 = [ver_opts.yticks, ver_opts.xticks], ver_opts.xticks = ref30[0], ver_opts.yticks = ref30[1];
    ref31 = [ver_opts.nyticks, ver_opts.nxticks], ver_opts.nxticks = ref31[0], ver_opts.nyticks = ref31[1];
  }
  verslice = d3panels.panelframe(ver_opts);
  g_heatmap = svg.append("g").attr("id", "heatmap");
  myheatmap(g_heatmap, data);
  g_horslice = svg.append("g").attr("id", "horslice").attr("transform", "translate(0," + htop + ")");
  horslice(g_horslice);
  g_verslice = svg.append("g").attr("id", "verslice").attr("transform", "translate(" + wleft + ",0)");
  verslice(g_verslice);
  formatX = d3panels.formatAxis(data.x);
  formatY = d3panels.formatAxis(data.y);
  cells = myheatmap.cells().on("mouseover", function(d, i) {
    g_verslice.select("g.title text").text("X = " + (formatX(d.x)));
    g_horslice.select("g.title text").text("Y = " + (formatY(d.y)));
    plotVer(d.xindex);
    return plotHor(d.yindex);
  }).on("mouseout", function(d, i) {
    g_verslice.select("g.title text").text("");
    return g_horslice.select("g.title text").text("");
  });
  vercurve = null;
  horcurve = null;
  plotHor = function(j) {
    if (horcurve != null) {
      horcurve.remove();
    }
    horcurve = d3panels.add_curves({
      linecolor: linecolor,
      linewidth: linewidth
    });
    return horcurve(horslice, {
      x: [data.x],
      y: [z_transpose[j]]
    });
  };
  plotVer = function(i) {
    if (vercurve != null) {
      vercurve.remove();
    }
    vercurve = d3panels.add_curves({
      linecolor: linecolor,
      linewidth: linewidth
    });
    if (flip_vert_slice) {
      return vercurve(verslice, {
        y: [data.y],
        x: [data.z[i]]
      });
    } else {
      return vercurve(verslice, {
        x: [data.y],
        y: [data.z[i]]
      });
    }
  };
  if (chartOpts.caption != null) {
    return d3.select(widgetdiv).insert("p").attr("class", "caption").text(chartOpts.caption);
  }
};
