// Generated by CoffeeScript 1.12.6
var add_search_box, iplotMap;

iplotMap = function(widgetdiv, data, chartOpts) {
  var axispos, chartdivid, clean_marker_name, div, height, horizontal, linecolor, linecolorhilit, linewidth, margin, markerSelect, martip, mychart, nyticks, rectcolor, ref, ref1, ref10, ref11, ref12, ref13, ref14, ref15, ref16, ref17, ref18, ref19, ref2, ref3, ref4, ref5, ref6, ref7, ref8, ref9, selectedMarker, shiftStart, svg, tickwidth, title, titlepos, widgetdivid, width, xlab, xlineOpts, ylab, ylim, yticks;
  width = (ref = chartOpts != null ? chartOpts.width : void 0) != null ? ref : 1000;
  height = (ref1 = chartOpts != null ? chartOpts.height : void 0) != null ? ref1 : 600;
  margin = (ref2 = chartOpts != null ? chartOpts.margin : void 0) != null ? ref2 : {
    left: 60,
    top: 40,
    right: 100,
    bottom: 40,
    inner: 10
  };
  axispos = (ref3 = chartOpts != null ? chartOpts.axispos : void 0) != null ? ref3 : {
    xtitle: 25,
    ytitle: 30,
    xlabel: 5,
    ylabel: 5
  };
  titlepos = (ref4 = chartOpts != null ? chartOpts.titlepos : void 0) != null ? ref4 : 20;
  ylim = (ref5 = chartOpts != null ? chartOpts.ylim : void 0) != null ? ref5 : null;
  nyticks = (ref6 = chartOpts != null ? chartOpts.nyticks : void 0) != null ? ref6 : 5;
  yticks = (ref7 = chartOpts != null ? chartOpts.yticks : void 0) != null ? ref7 : null;
  xlineOpts = (ref8 = chartOpts != null ? chartOpts.xlineOpts : void 0) != null ? ref8 : {
    color: "#cdcdcd",
    width: 5
  };
  tickwidth = (ref9 = chartOpts != null ? chartOpts.tickwidth : void 0) != null ? ref9 : 10;
  rectcolor = (ref10 = chartOpts != null ? chartOpts.rectcolor : void 0) != null ? ref10 : "#E6E6E6";
  linecolor = (ref11 = chartOpts != null ? chartOpts.linecolor : void 0) != null ? ref11 : "slateblue";
  linecolorhilit = (ref12 = chartOpts != null ? chartOpts.linecolorhilit : void 0) != null ? ref12 : "Orchid";
  linewidth = (ref13 = chartOpts != null ? chartOpts.linewidth : void 0) != null ? ref13 : 3;
  title = (ref14 = chartOpts != null ? chartOpts.title : void 0) != null ? ref14 : "";
  xlab = (ref15 = chartOpts != null ? chartOpts.xlab : void 0) != null ? ref15 : "Chromosome";
  ylab = (ref16 = chartOpts != null ? chartOpts.ylab : void 0) != null ? ref16 : "Position (cM)";
  shiftStart = (ref17 = chartOpts != null ? chartOpts.shiftStart : void 0) != null ? ref17 : false;
  horizontal = (ref18 = chartOpts != null ? chartOpts.horizontal : void 0) != null ? ref18 : false;
  chartdivid = (ref19 = chartOpts != null ? chartOpts.chartdivid : void 0) != null ? ref19 : 'chart';
  widgetdivid = d3.select(widgetdiv).attr('id');
  margin = d3panels.check_listarg_v_default(margin, {
    left: 60,
    top: 40,
    right: 100,
    bottom: 40,
    inner: 10
  });
  axispos = d3panels.check_listarg_v_default(axispos, {
    xtitle: 25,
    ytitle: 30,
    xlabel: 5,
    ylabel: 5
  });
  mychart = d3panels.mapchart({
    height: height,
    width: width,
    margin: margin,
    axispos: axispos,
    titlepos: titlepos,
    ylim: ylim,
    yticks: yticks,
    nyticks: nyticks,
    xlineOpts: xlineOpts,
    tickwidth: tickwidth,
    rectcolor: rectcolor,
    linecolor: linecolor,
    linecolorhilit: linecolorhilit,
    linewidth: linewidth,
    title: title,
    xlab: xlab,
    ylab: ylab,
    horizontal: horizontal,
    shiftStart: shiftStart,
    tipclass: widgetdivid
  });
  div = d3.select(widgetdiv);
  mychart(div.select("svg"), data);
  svg = mychart.svg();
  martip = d3.tip().attr('class', "d3-tip " + widgetdivid).html(function(d) {
    var pos;
    pos = d3.format(".1f")(data.pos[data.marker.indexOf(d)]);
    return d + " (" + pos + ")";
  }).direction(function() {
    if (horizontal) {
      return 'n';
    }
    return 'e';
  }).offset(function() {
    if (horizontal) {
      return [-10, 0];
    }
    return [0, 10];
  });
  svg.call(martip);
  clean_marker_name = function(markername) {
    return markername.replace(".", "\\.").replace("#", "\\#").replace("/", "\\/");
  };
  selectedMarker = "";
  $("div#markerinput_" + widgetdivid).submit(function() {
    var line, newSelection;
    newSelection = document.getElementById("marker_" + widgetdivid).value;
    event.preventDefault();
    if (selectedMarker !== "") {
      div.select("line#" + (clean_marker_name(selectedMarker))).attr("stroke", linecolor);
      martip.hide();
    }
    if (newSelection !== "") {
      if (data.marker.indexOf(newSelection) >= 0) {
        selectedMarker = newSelection;
        line = div.select("line#" + (clean_marker_name(selectedMarker))).attr("stroke", linecolorhilit);
        martip.show(line.datum(), line.node());
        div.select("a#currentmarker").text("");
        return true;
      } else {
        div.select("a#currentmarker").text("Marker \"" + newSelection + "\" not found");
      }
    }
    return false;
  });
  $("input#marker_" + widgetdivid).autocomplete({
    autoFocus: true,
    source: function(request, response) {
      var matches;
      matches = $.map(data.marker, function(tag) {
        if (tag.toUpperCase().indexOf(request.term.toUpperCase()) === 0) {
          return tag;
        }
      });
      return response(matches);
    },
    select: function(event, ui) {
      $("input#marker_" + widgetdivid).val(ui.item.label);
      return $("input#submit_" + widgetdivid).submit();
    }
  });
  $("input#marker_" + widgetdivid).each(function() {
    $("div.searchbox#markerinput_" + widgetdivid).addClass('inactive');
    return $(this).data('default', $(this).val()).focus(function() {
      $("div.searchbox#markerinput_" + widgetdivid).removeClass('inactive');
      if ($(this).val() === $(this).data('default') || $(this).val() === '') {
        return $(this).val('');
      }
    }).blur(function() {
      if ($(this).val() === '') {
        $("div.searchbox#markerinput_" + widgetdivid).addClass('inactive');
        return $(this).val($(this).data('default'));
      }
    });
  });
  markerSelect = mychart.markerSelect();
  markerSelect.on("mouseover", function(d) {
    if (selectedMarker !== "") {
      if (selectedMarker !== d) {
        div.select("line#" + (clean_marker_name(selectedMarker))).attr("stroke", linecolor);
      }
      return martip.hide();
    }
  });
  if (chartOpts.caption != null) {
    return d3.select(widgetdiv).insert("p").attr("class", "caption").text(chartOpts.caption);
  }
};

add_search_box = function(widgetdiv) {
  var div, form, widgetdivid;
  div = d3.select(widgetdiv);
  widgetdivid = div.attr("id");
  form = div.append("div").attr("class", "searchbox").attr("id", "markerinput_" + widgetdivid).append("form").attr("name", "markerinput_" + widgetdivid);
  form.append("input").attr("id", "marker_" + widgetdivid).attr("type", "text").attr("value", "Marker name").attr("name", "marker");
  form.append("input").attr("type", "submit").attr("id", "submit_" + widgetdivid).attr("value", "Submit");
  return form.append("a").attr("id", "currentmarker");
};
