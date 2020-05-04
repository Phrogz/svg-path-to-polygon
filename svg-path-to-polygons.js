module.exports = {
  pathDataToPolys:svgPathToPolygons,
  compare:compare
};
const { parseSVG, makeAbsolute } = require('svg-path-parser');

function svgPathToPolygons(svgPathString, opts) {
  return new SvgPathToPolygons(svgPathString, opts);
}

function SvgPathToPolygons(svgPathString, opts={}) {
  if (!opts.tolerance) opts.tolerance=1;
  this.opts = opts;
  this.tolerance2 = opts.tolerance*opts.tolerance;
  this.polys = [];
  this.poly = [];

  makeAbsolute(parseSVG(svgPathString)).forEach(this._interpretPath.bind(this));

  return this.polys;
}

// convert path to points
SvgPathToPolygons.prototype._interpretPath = function (cmd) {
  switch(cmd.code) {
    case 'M':
      this.polys.push(this.poly=[]);
      // intentional flow-through
    case 'L':
    case 'H':
    case 'V':
    case 'Z':
      this._add(cmd.x,cmd.y);
      if (cmd.code==='Z')
        this.poly.closed = true;
    break;
    case 'A':
      ellipticArcToCubicBezierCurves(
        [cmd.x0, cmd.y0],
        [cmd.x, cmd.y],
        [cmd.rx, cmd.ry],
        cmd.xAxisRotation,
        cmd.largeArc,
        cmd.sweep).forEach((curve, i) => {
        const cp0 = [curve[0][0], curve[0][1]];
        const cp1 = [curve[2][0], curve[2][1]];
        const cp2 = [curve[3][0], curve[3][1]];
        const cp3 = [curve[1][0], curve[1][1]];

        if (i==0)
          this._add(cp0[0],cp0[1]);

        this._sampleCubicBezier(cp0[0], cp0[1],
                                cp1[0], cp1[1],
                                cp2[0], cp2[1],
                                cp3[0], cp3[1]);
        this._add(cp3[0],cp3[1]);
      });

      this._add(cmd.x,cmd.y);
      break;
    case 'C':
      this._sampleCubicBezier(cmd.x0,cmd.y0,cmd.x1,cmd.y1,cmd.x2,cmd.y2,cmd.x,cmd.y);

      this._add(cmd.x,cmd.y);
    break;

    case 'Q':
      const [cp0, cp1, cp2, cp3] = quadraticToCubicBezier(cmd);

      this._sampleCubicBezier(cp0[0], cp0[1],
                              cp1[0], cp1[1],
                              cp2[0], cp2[1],
                              cp3[0], cp3[1]);
      this._add(cp3[0],cp3[1]);
    break;

    case 'S':
      let x1=0, y1=0;
      if (this.prev) {
        if (this.prev.code==='C') {
          x1 = this.prev.x*2 - this.prev.x2;
          y1 = this.prev.y*2 - this.prev.y2;
        } else {
          x1 = this.prev.x;
          y1 = this.prev.y;
        }
      }
      this._sampleCubicBezier(cmd.x0,cmd.y0,x1,y1,cmd.x2,cmd.y2,cmd.x,cmd.y);
      this._add(cmd.x,cmd.y);
    break;

    default:
      console.error('Our deepest apologies, but '+cmd.command+' commands ('+cmd.code+') are not yet supported.');
      process.exit(2);
  }

  this.prev = cmd;
};

// args
// start point     : [x0, y0]
// control point 1 : [x1, y1]
// control point 2 : [x2, y2]
// final point     : [x3, y3]
// http://antigrain.com/research/adaptive_bezier/
SvgPathToPolygons.prototype._sampleCubicBezier = function (x0, y0, x1, y1, x2, y2, x3, y3) {
  // Calculate all the mid-points of the line segments
  const x01   = (x0 + x1) / 2,
        y01   = (y0 + y1) / 2,
        x12   = (x1 + x2) / 2,
        y12   = (y1 + y2) / 2,
        x23   = (x2 + x3) / 2,
        y23   = (y2 + y3) / 2,
        x012  = (x01 + x12) / 2,
        y012  = (y01 + y12) / 2,
        x123  = (x12 + x23) / 2,
        y123  = (y12 + y23) / 2,
        x0123 = (x012 + x123) / 2,
        y0123 = (y012 + y123) / 2;

  // Try to approximate the full cubic curve by a single straight line
  const dx = x3-x0,
        dy = y3-y0;

  const d1 = Math.abs(((x1-x3)*dy - (y1-y3)*dx)),
        d2 = Math.abs(((x2-x3)*dy - (y2-y3)*dx));

  if (((d1+d2)*(d1+d2)) < (this.tolerance2 * (dx*dx + dy*dy))) {
    this._add(x0123,y0123);
    return ;
  } else { // Continue subdivision
    this._sampleCubicBezier(x0, y0, x01, y01, x012, y012, x0123, y0123);
    this._sampleCubicBezier(x0123, y0123, x123, y123, x23, y23, x3, y3);
  }
};

SvgPathToPolygons.prototype._add = function (x,y) {
  const decimals = this.opts.decimals;
  if (decimals && decimals>=0) {
    x = x.toFixed(decimals)*1;
    y = y.toFixed(decimals)*1;
  }
  this.poly.push([x,y]);
};

//
// helping function to convert some path to cubic bezier curves
//

// convert quadratic cuve to cubic bezier one
// https://stackoverflow.com/questions/3162645/convert-a-quadratic-bezier-to-a-cubic-one
//
// arguments cmd
//
// return
// start point     : [x0, y0]
// control point 1 : [x1, y1]
// control point 2 : [x2, y2]
// final point     : [x3, y3]
function quadraticToCubicBezier(cmd) {
  let qp0 = [cmd.x0, cmd.y0];
  let qp1 = [cmd.x1, cmd.y1];
  let qp2 = [cmd.x, cmd.y];

  let cp0 = [qp0[0], qp0[1]];
  let cp3 = [qp2[0], qp2[1]];
  let cp1 = [
    qp0[0] + 2/3. * (qp1[0] - qp0[0]),
    qp0[1] + 2/3. * (qp1[1] - qp0[1])
  ];
  let cp2 = [
    qp2[0] + 2/3. * (qp1[0] - qp2[0]),
    qp2[1] + 2/3. * (qp1[1] - qp2[1])
  ];

  return [cp0, cp1, cp2, cp3];
}

// convert and ellipticArc to multiple cubic bezier curve
//
// source
// https://mortoray.com/2017/02/16/rendering-an-svg-elliptical-arc-as-bezier-curves/
// https://github.com/fuse-open/fuselibs/blob/master/Source/Fuse.Drawing.Surface/SurfaceUtil.uno
//
// arguments
// p1 = start point [x, y]
// p2 = final point [x, y]
// r_ = [rx, ry]
// xAngle = ...
// flagA = ...
// flagS = ...
//
// will
// return [ [
//   [start point],
//   [final point],
//   [control point 1],
//   [control point 2],
// ], ...]
var ellipticArcToCubicBezierCurves;

(function () {
  const _zeroTolerance = 1e-05;

  ellipticArcToCubicBezierCurves = function (p1, p2, r_, xAngle, flagA, flagS) {
    const curveInfo = endpointToCenterArcParams(p1, p2, r_, xAngle, flagA, flagS);
    const newR = curveInfo[0];
    const center = curveInfo[1];
    const angle = curveInfo[2];

    return ellipticArcToBezierCurve(center, newR, xAngle, angle[0], angle[1]);
  }

  function ellipticArcToBezierCurve(center, radius, xAngle, startAngle, deltaAngle) {
    const curves = [];
    var s = startAngle;
    const e = s + deltaAngle;
    const neg = e < s;
    const sign = neg ? -1 : 1;
    var remain = Math.abs(e - s);

    var prev = ellipticArcPoint(center, radius, xAngle, s);

    while( remain > _zeroTolerance ) {
      const step = Math.min( remain, Math.PI / /* 4 */ 8 );
      const signStep = step * sign;


      const p1 = prev;
      const p2 = ellipticArcPoint( center, radius, xAngle, s + signStep );

      const alphaT = Math.tan(signStep / 2);
      const alpha = Math.sin(signStep) * (Math.sqrt(4 + 3 * alphaT * alphaT)- 1) / 3;

      const el1 = ellipticArcDerivative(center, radius, xAngle, s);
      const el2 = ellipticArcDerivative(center, radius, xAngle, s + signStep);

      const q1x = p1[0] + alpha * el1[0];
      const q2x = p2[0] - alpha * el2[0];
      const q1y = p1[1] + alpha * el1[1];
      const q2y = p2[1] - alpha * el2[1];

      const q1 = [q1x, q1y];
      const q2 = [q2x, q2y];

      curves.push([p1, p2, q1, q2]);

      s += signStep;
      remain -= step;
      prev = p2;
    }
    return curves;
  }


  /*
  Equations from:
    Drawing an elliptical arc using polylines, quadratic or cubic Bézier curves
      by L. Maisonobe
  http://www.spaceroots.org/documents/ellipse/elliptical-arc.pdf
  */
  function ellipticArcPoint(c, r, xAngle,  t) {
    return [
      c[0] + r[0] * Math.cos(xAngle) * Math.cos(t) - r[1] * Math.sin(xAngle) * Math.sin(t),
      c[1] + r[0] * Math.sin(xAngle) * Math.cos(t) + r[1] * Math.cos(xAngle) * Math.sin(t)
    ];
  }

  function ellipticArcDerivative(c, r, xAngle, t) {
    return [
      -r[0] * Math.cos(xAngle) * Math.sin(t) - r[1] * Math.sin(xAngle) * Math.cos(t),
      -r[0] * Math.sin(xAngle) * Math.sin(t) + r[1] * Math.cos(xAngle) * Math.cos(t)
    ];
  }

  function Math_fmod(a,b) { return Number((a - (Math.floor(a / b) * b))); };
  function Math_clamp(val, min, max) { return Math.min(Math.max(min, val), max); }

  /**
    Perform the endpoint to center arc parameter conversion as detailed in the SVG 1.1 spec.
    F.6.5 Conversion from endpoint to center parameterization

    @param r must be a ref in case it needs to be scaled up, as per the SVG spec

    arguments
    p1 = start point [x, y]
    p2 = final point [x, y]
    r_ = [rx, ry]

    return [r_[2],  center[2], angles[2]]
    */
  function endpointToCenterArcParams(p1, p2, r_, xAngle, flagA, flagS) {
    let rX = Math.abs(r_[0]);
    let rY = Math.abs(r_[0]);

    //(F.6.5.1)
    let dx2 = (p1[0] - p2[0]) / 2.0;
    let dy2 = (p1[1] - p2[1]) / 2.0;
    let x1p = Math.cos(xAngle)*dx2 + Math.sin(xAngle)*dy2;
    let y1p = -Math.sin(xAngle)*dx2 + Math.cos(xAngle)*dy2;

    //(F.6.5.2)
    let rxs = rX * rX;
    let rys = rY * rY;
    let x1ps = x1p * x1p;
    let y1ps = y1p * y1p;
    // check if the radius is too small `pq < 0`, when `dq > rxs * rys` (see below)
    // cr is the ratio (dq : rxs * rys)
    let cr = x1ps/rxs + y1ps/rys;
    if (cr > 1) {
      //scale up rX,rY equally so cr == 1
      var s = Math.sqrt(cr);
      rX = s * rX;
      rY = s * rY;
      rxs = rX * rX;
      rys = rY * rY;
    }
    let dq = (rxs * y1ps + rys * x1ps);
    let pq = (rxs*rys - dq) / dq;
    let q = Math.sqrt( Math.max(0,pq) ); //use Max to account for float precision
    if (flagA === flagS)
      q = -q;
    let cxp = q * rX * y1p / rY;
    let cyp = - q * rY * x1p / rX;

    //(F.6.5.3)
    let cx = Math.cos(xAngle)*cxp - Math.sin(xAngle)*cyp + (p1[0] + p2[0])/2;
    let cy = Math.sin(xAngle)*cxp + Math.cos(xAngle)*cyp + (p1[1] + p2[1])/2;

    //(F.6.5.5)
    let theta = svgAngle( 1,0, (x1p-cxp) / rX, (y1p - cyp)/rY );
    //(F.6.5.6)
    let delta = svgAngle(
      (x1p - cxp)/rX, (y1p - cyp)/rY,
      (-x1p - cxp)/rX, (-y1p-cyp)/rY);

    delta = Math_fmod(delta, Math.PI * 2 );

    if (!flagS)
      delta -= 2 * Math.PI;

    r_ = [rX, rY];
    let c = [cx, cy];
    let angles = [theta, delta];

    return [r_, c, angles];
  }

  function _2d_vector_dot(a, b) {
    return a[0]*b[0] + a[1]*b[1];
  }

  function _2d_vector_length(a) {
    return Math.sqrt(_2d_vector_dot(a, a));
  }

  function svgAngle(ux, uy, vx, vy) {
    var u = [ux, uy];
    var v = [vx, vy];
    //(F.6.5.4)
    var dot = _2d_vector_dot(u, v);
    var len = _2d_vector_length(u) * _2d_vector_length(v);
    var ang = Math.acos( Math_clamp(dot / len,-1,1) ); //floating point precision, slightly over values appear
    if ( (u[0]*v[1] - u[1]*v[0]) < 0)
      ang = -ang;
    return ang;
  }
})();

// OMG YOU FOUND THE SECRET UNDOCUMENTED FEATURE
function compare(pathData,opts={}) {
  var polys = svgPathToPolygons(pathData,opts);
  var minX=Infinity, maxX=-Infinity, minY=Infinity, maxY=-Infinity;
  polys.forEach(poly => {
    poly.forEach(pt => {
      if (pt[0]<minX) minX=pt[0];
      if (pt[1]<minY) minY=pt[1];
      if (pt[0]>maxX) maxX=pt[0];
      if (pt[1]>maxY) maxY=pt[1];
    });
  });
  let dx=maxX-minX, dy=maxY-minY;
  console.log(`
<svg xmlns="http://www.w3.org/2000/svg" width="${dx}px" height="${dy}px" viewBox="${minX} ${minY} ${dx*2} ${dy}">
<style>path,polygon,polyline { fill-opacity:0.2; stroke:black }</style>
<path d="${pathData}"/>
<g transform="translate(${dx},0)">
${polys.map(poly => `  <${poly.closed ? 'polygon' : 'polyline'} points="${poly.join(' ')}"/>`).join("\n")}
</g>
</svg>
  `.trim());
};
