const { assert } = require('chai');
const { pathDataToPolys } = require('../svg-path-to-polygons.js');

const expectations = [
	{
		m:'linear open',
		d:'M5,7 10,20 30,40',
		o:[ [[5,7],[10,20],[30,40]] ],
		closed: [ false ]
	},
	{
		m:'linear closed',
		d:'M5,7 10,20 30,40 z',
		o:[ [[5,7],[10,20],[30,40],[5,7]] ],
		closed: [ true ]
	},
	{
		m:'linear two path',
		d:'M5,7 10,20 30,40 z M100,100',
		o:[
			[[5,7],[10,20],[30,40],[5,7]],
			[[100,100]]
		],
		closed: [ true, false ]
	},
	{
		m:'linear two path with relative',
		d:'M5,7 10,20 30,40 V10 H20 v-10 h-10z m5,2 l5,5 5,5',
		o:[
			[[5,7],[10,20],[30,40],[30,10],[20,10],[20,0],[10,0],[5,7]],
			[[10,9],[15,14],[20,19]]
		],
		closed: [ true, false ]
	},
	{
		m:'unsplit curve',
		d:'M5,15 c5.5,0 10-4.5 10,-10 h10',
		o:[ [[5,15],[12.0625,12.0625],[15,5],[25,5]] ],
		closed: [false],
		tolerance:1000
	},
	{
		m:'default curve split',
		d:'M5,15 c5.5,0 10-4.5 10,-10 h10',
		o:[[[5,15],[7,14.8],[10.6,13.3],[13.3,10.6],[14.8,7],[15,5],[25,5]]],
		closed: [false],
		decimals:1
	},
	{
		m:'degenerate curve',
		d:'M0,0 c0,0 0,0 0,0',
		o:[[[0,0],[0,0]]],
		closed: [false]
	},

];

expectations.forEach(ex => {
	let o = pathDataToPolys(ex.d, {tolerance:ex.tolerance, decimals:ex.decimals});
	assert.deepEqual(o, ex.o, ex.m+', path data');
	assert.deepEqual(o.map( poly=>!!poly.closed), ex.closed, ex.m+', properly closed' );
});

// TODO: handle parse errors
