/*

 Chromanin 
 
 Based on public domain pascal code of Texture v.0.5 (C)2K2 by CARSTEN WAECHTER http://ainc.de
 
 JavaScript version Copyright (C) 2012 Alejandro Mosquera <amsqr2@gmail.com>
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, version 3 of the License.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

function Chromanin(texsize,strcanvas)
{
this.MAX_LAYERS=5;
this.TEMPL =5;
this.layerSizeX=texsize;
this.layerSizeY=texsize;
this.andLayerSizeX=texsize-1;
this.andLayerSizeY=texsize-1;
this.layers = []; // array de canvas
this.imageData=null;
this.seedValue=0;
this.outcanvas = strcanvas;

}


function RGBA()
{
	this.r=0;
	this.g=0;
	this.b=0;
	this.a=255;
};

Chromanin.prototype.setPixel= function setPixel(x, r, g, b, a) {
    index = (x);
    this.imageData.data[index+0] = r;
    this.imageData.data[index+1] = g;
    this.imageData.data[index+2] = b;
    this.imageData.data[index+3] = a;
}

Chromanin.prototype.writeCanvas= function writeCanvas()
{
    var skel = "ch_canvas_tmp";
	
	for (var i=0; i<this.MAX_LAYERS; i++)
	{
	    ref = null;
		if (i==4) 
		{
			ref = document.getElementById(this.outcanvas);
		}
		else
		{
			ref = document.getElementById(skel+i);
			//alert(skel+i);
		}
		
		if (i<=4)
		{
		
		ctx = ref.getContext('2d');
		this.imageData = ctx.createImageData(this.layerSizeX, this.layerSizeY);
		//alert(i);
		for (var x=0; x<(this.layerSizeX * this.layerSizeY); x+=1)
		{
		    index=4*x;
		    this.imageData.data[index] =  Math.round(this.layers[i][x].r);
			this.imageData.data[index+1] =  Math.round(this.layers[i][x].g);
			this.imageData.data[index+2] =  Math.round(this.layers[i][x].b);
			this.imageData.data[index+3] = Math.round(this.layers[i][x].a)
			//setPixel(x, Math.round(this.layers[i][x].r),Math.round(this.layers[i][x].g),Math.round(this.layers[i][x].b),Math.round(this.layers[i][x].a));
		}
		
		
		
		
		
	    ctx.putImageData(this.imageData, 0, 0); 
	  }
	
	}
    
}





Chromanin.prototype.CreateTexture= function CreateTexture()
{
	


}



Chromanin.prototype.initlayers= function initlayers(sizeX, sizeY)
{

	/*
	for (var i=0; i<this.MAX_LAYERS+1; i++)
	{
		var skel = "ch_canvas_tmp";
		
	    ref = document.createElement('canvas');
		ref.id = skel+i;
		ref.width = sizeX; 
        ref.height =sizeY; 
		document.body.appendChild(ref);
		
		ref = document.getElementById(skel+i);
		//alert(ref);
	}
	
	*/
	
	for (var i=0; i<this.MAX_LAYERS+1; i++)
	{
	
		this.layers[i]=[];
		for (var x=0; x<(sizeX * sizeY); x+=1)
		{
			this.layers[i][x] = new RGBA();
		}
	
	}
	
	
	
	this.layerSizeX=sizeX;
	this.layerSizeY=sizeY;
	this.andLayerSizeX=sizeX-1;
	this.andLayerSizeY=sizeY-1;
}


Chromanin.prototype.deInitLayers= function deInitLayers()
{
	
	
}


Chromanin.prototype.getLayer= function getLayer(l)
{
	return this.layers[l];
}


function myRandom()
{
        this.seedValue = (this.seedValue * 22695477 + 1) & 0xffffff;
        return (this.seedValue >> 16) & 0x7fff;
        //return Math.floor(Math.random()*301);
}


function cosineInterpolate(v, x, y)
{
	var f1, f2, mf1, mf2, g0, g1, g2, g3;
	color= new RGBA();
	mf1=(1-Math.cos(x*Math.PI))*0.5;
	mf2=(1-Math.cos(y*Math.PI))*0.5;
	f1=1-mf1;
	f2=1-mf2;
	g0=f1*f2;
	g1=mf1*f2;
	g2=f1*mf2;
	g3=mf1*mf2;

	color.r=(v[0].r*g0+v[1].r*g1+v[2].r*g2+v[3].r*g3);
	color.g=(v[0].g*g0+v[1].g*g1+v[2].g*g2+v[3].g*g3);
	color.b=(v[0].b*g0+v[1].b*g1+v[2].b*g2+v[3].b*g3);

	return color;
}


function mySqrtInt(n)
{
	var root, tryr;

	if (n>=65025)		// 255*255
		return 255;

	root=0;
	for (var i=15; i>=0; i--)
	{
		tryr=root+(1<<i);
		if (n>(tryr<<i))
		{
			n-=tryr<<i;
			root=root|(2<<i);
		}
	}
	return root>>1;
}

function fmod(a, b)
{
	var x = Math.floor(a/b);
	return a - b*x;
}




function hsv2rgb( h,  s, v, r, g, b)
{
	var f, p, q, t;
	var i;

	while (h<0) h+=360;
	while (h>=360) h-=360;

	if (s==0)
	{
		r=v;
		g=v;
		b=v;
		}
	else {
		h/=60;
		i=getlo(h);
		f=fmod(h, 1);
		f=(h-i);
		p=v*(1-s);
		q=v*(1-(s*f));
		t=v*(1-(s*(1-f)));
		switch(i)
		{
			case 0:	r=v; g=t; b=p; break;
			case 1:	r=q; g=v; b=p; break;
			case 2:	r=p; g=v; b=t; break;
			case 3:	r=p; g=q; b=v; break;
			case 4:	r=t; g=p; b=v; break;
			case 5:	r=v; g=p; b=q; break;
		}
	}
}


function rgb2hsv(r, g, b, h, s, v)
{
	var maxR, maxG, maxB, delta;
	var mx=Math.max(r, Math.max(g, b));
	var mn=Math.min(r, Math.min(g, b));

	v=mx;
	s=0;

	if (mx!=0)
		s=(mx-mn)/mx;
	if (s==0)
		h=-1;
	else {
		delta=mx-mn;
		maxR=mx-r;
		maxG=mx-g;
		maxB=mx-b;
		if (r==mx)
			h=(maxB-maxG)/delta;
		else
			if (g==mx)
				h=2+(maxR-maxB)/delta;
			else
				h=4+(maxG-maxR)/delta;
		h*=60;
		while (h<0) h+=360;
		while (h>=360) h-=360;
	}
}


Chromanin.prototype.getBilerpixel= function getBilerpixel(l, x, y)
{
	var corner = [];
	 corner[0] =new RGBA();
	 corner[1] =new RGBA();
	 corner[2] =new RGBA();
	 corner[3] =new RGBA();
	var xi, yi, xip, xip1;
	var yip, yip1;

	xi=x;
	yi=y;
	xip=xi&this.andLayerSizeX;
	xip1=(xi+1)&this.andLayerSizeX;
	yip1=((yi+1)&this.andLayerSizeY)*this.layerSizeX;
	yip=(yi&this.andLayerSizeY)*this.layerSizeX;

	corner[0]=this.layers[l][xip+yip];
	corner[1]=this.layers[l][xip1+yip];
	corner[2]=this.layers[l][xip+yip1];
	corner[3]=this.layers[l][xip1+yip1];

	return cosineInterpolate(corner, (x-xi), (y-yi));
}

Chromanin.prototype.subPlasma= function subPlasma(l, dist, seed, amplitude,  rgb)
{
	var x, y;
	var offset, offset2;
	var corner = [];
	 corner[0] =new RGBA();
	 corner[1] =new RGBA();
	 corner[2] =new RGBA();
	 corner[3] =new RGBA();
	var oodist;

	if (seed!=0)
		this.seedValue=seed;

	for (y=0; y<this.layerSizeY; y+=dist)
		for (x=0; x<this.layerSizeX; x+=dist)
		{
			offset=y*this.layerSizeX+x;
			this.layers[l][offset].r=this.layers[l][offset].g=this.layers[l][offset].b=getlo(myRandom())&(amplitude-1);
			if (rgb)
			{
				this.layers[l][offset].g=getlo(myRandom())&(amplitude-1);
				this.layers[l][offset].b=getlo(myRandom())&(amplitude-1);
			}
		}

	if (dist<1)
		return;

	oodist=1/dist;

	for (y=0; y<this.layerSizeY; y+=dist)
	{
		offset=y*this.layerSizeX;
		offset2=((y+dist)&this.andLayerSizeY)*this.layerSizeX;
		for (x=0; x<this.layerSizeX; x+=dist)
		{
			corner[0]=this.layers[l][x+offset];
			corner[1]=this.layers[l][((x+dist)&this.andLayerSizeX)+offset];
			corner[2]=this.layers[l][x+offset2];
			corner[3]=this.layers[l][((x+dist)&this.andLayerSizeX)+offset2];
			for (var b=0; b<dist; b++)
				for (var a=0; a<dist; a++)
					this.layers[l][x+a+(y+b)*this.layerSizeX]=cosineInterpolate(corner, oodist*a, oodist*b);
		}
	}
}


Chromanin.prototype.sinePlasma= function sinePlasma(l, dx,  dy,  amplitude)
{
	amplitude/=256;
	for (var y=0; y<this.layerSizeY; y++)
		for (var x=0; x<this.layerSizeX; x++)
			this.layers[l][x+y*this.layerSizeX].r=this.layers[l][x+y*this.layerSizeX].g=this.layers[l][x+y*this.layerSizeX].b=((63.5*Math.sin(dx*x)+127+63.5*Math.sin(dy*y))*amplitude);
	
}


Chromanin.prototype.perlinNoise= function perlinNoise( l,  dist,  seed,  amplitude,  persistence,  octaves,  rgb)
{
	 var r;

	this.subPlasma(l, dist, seed, 1, rgb);
	for ( var i=0; i<octaves-2; i++)
	{
		amplitude=(amplitude*persistence)>>8;
		if (amplitude<=0) break;
		dist=dist>>1;
		if (dist<=0) break;
		this.subPlasma(this.TEMPL, dist, 0, amplitude, rgb);
		for (var v=0; v<this.layerSizeX*this.layerSizeY; v++)
		{
			r=this.layers[l][v].r+this.layers[this.TEMPL][v].r;
			if (r>255)
				this.layers[l][v].r=255;
			else
				this.layers[l][v].r=r;
			r=this.layers[l][v].g+this.layers[this.TEMPL][v].g;
			if (r>255)
				this.layers[l][v].g=255;
			else
				this.layers[l][v].g=r;
			r=this.layers[l][v].b+this.layers[this.TEMPL][v].b;
			if (r>255)
				this.layers[l][v].b=255;
			else
				this.layers[l][v].b=r;
		}
	}
}


Chromanin.prototype.particle= function particle( l,  f)
{
	 var r, nx, ny;
	 var offset;
	f*=180;
	for (var y=0; y<this.layerSizeY; y++)
	{
		ny=y/(this.layerSizeY>>1)-1;
		for (var x=0; x<this.layerSizeX; x++)
		{
			offset=y*this.layerSizeX+x;
			nx=x/(this.layerSizeX>>1)-1;
			r=255-f*Math.sqrt(nx*nx+ny*ny);
			if (r<0) r=0;
			if (r>255) r=255;
			this.layers[l][offset].r=this.layers[l][offset].g=this.layers[l][offset].b=(r);
		}
	}
}


Chromanin.prototype.colorLayer= function colorLayer( l,  r,  g,  b)
{
	var color= new RGBA();
	color.r=r; color.g=g; color.b=b;
	for (var v=0; v<this.layerSizeX*this.layerSizeY; v++)
		this.layers[l][v]=color;
}


Chromanin.prototype.checkerBoardLayer= function checkerBoardLayer( l,  dx,  dy,  r1,  g1,  b1,  r2,  g2,  b2)
{
	var col1=new RGBA();
	var col2=new RGBA();
	col1.r=r1; col1.g=g1; col1.b=b1;
	col2.r=r2; col2.g=g2; col2.b=b2;
	for (var y=0; y<this.layerSizeY; y++)
		for (var x=0; x<this.layerSizeX; x++)
			if (((y/dy)&1)^((x/dx)&1))
				this.layers[l][y*this.layerSizeX+x]=col1;
			else
				this.layers[l][y*this.layerSizeX+x]=col2;
}


Chromanin.prototype.blobsLayer= function blobsLayer( l,  seed,  amount,  rgb)
{
	 var blobX =new Array(16);
 	 var blobY =new Array(16);
	 var blobR=new Array(16);
	 var blobG =new Array(16);
	 var blobB=new Array(16);
	 var r, g, b, sd, d, oosize;
	 var rx, ry;
	 var offset;

	 this.seedValue=seed;

	for ( v=0; v<amount; v++)
	{
		blobX[v]=getlo(getlo(myRandom())&this.andLayerSizeX);
		blobY[v]=getlo(getlo(myRandom())&this.andLayerSizeY);
		blobR[v]=(getlo(myRandom()&255)/255+0.1)*150;
		if (rgb==true)
		{
			blobG[v]=getlo((getlo(getlo(myRandom())&255)/255+0.1)*150);
			blobB[v]=getlo((getlo(getlo(myRandom())&255)/255+0.1)*150);
		}
	}

	oosize=3/(this.layerSizeX*this.layerSizeY);
	for (var y=0; y<this.layerSizeY; y++)
		for (var x=0; x<this.layerSizeX; x++)
		{
			offset=y*this.layerSizeX+x;
			r=g=b=0;
			for (v=0; v<amount; v++)
			{
				rx=blobX[v]-x;
				ry=blobY[v]-y;
				d=oosize*(rx*rx+ry*ry);
				sd=d*d;
				d=-0.444444*sd*d+1.888888*sd-2.444444*d+1;
				r+=d*blobR[v];
				g+=d*blobG[v];		// needn't be calculated if not rgb, but we do it for memory optimization (spares one if statement)
				b+=d*blobB[v];		// needn't be calculated if not rgb, but we do it for memory optimization (spares one if statement)
			}

			if (r<0) r=0;
			if (r>255) r=255;
			this.layers[l][offset].r=this.layers[l][offset].g=this.layers[l][offset].b=getlo(r);
			if (rgb)
			{
				if (g<0) g=0;
				if (g>255) g=255;
				this.layers[l][offset].g=getlo(g);
				if (b<0) b=0;
				if (b>255) b=255;
				this.layers[l][offset].b=getlo(b);
			}
		}
}


Chromanin.prototype.scaleLayerRGB= function scaleLayerRGB( src,  dest,  r,  g,  b)
{
	 var tr, tg, tb;

	for (var v=0; v<this.layerSizeX*this.layerSizeY; v++)
	{
		tr=this.layers[src][v].r*r;
		tg=this.layers[src][v].g*g;
		tb=this.layers[src][v].b*b;

		if (tr>=255) this.layers[dest][v].r=255; else if (tr<=0) this.layers[dest][v].r=0; else this.layers[dest][v].r=tr;
		if (tg>=255) this.layers[dest][v].g=255; else if (tg<=0) this.layers[dest][v].g=0; else this.layers[dest][v].g=tg;
		if (tb>=255) this.layers[dest][v].b=255; else if (tb<=0) this.layers[dest][v].b=0; else this.layers[dest][v].b=tb;
	}
}


Chromanin.prototype.scaleLayerHSV= function scaleLayerHSV( src,  dest,  h,  s,  v)
{
	 var th, ts, tv, tr, tg, tb;

	for (var k=0; k<this.layerSizeX*this.layerSizeY; k++)
	{
		rgb2hsv(this.layers[src][k].r, this.layers[src][k].g, this.layers[src][k].b, th, ts, tv);
		th*=h;
		ts*=s;
		tv*=v;
		if (ts>1) ts=1; else if (ts<0) ts=0;
		if (tv>255) tv=255; else if (tv<0) tv=0;
		hsv2rgb(th, ts, tv, tr, tg, tb);
		this.layers[dest][k].r=getlo(tr);
		this.layers[dest][k].g=getlo(tg);
		this.layers[dest][k].b=getlo(tb);
	}
}


Chromanin.prototype.adjustLayerRGB= function adjustLayerRGB( src,  dest,  r,  g,  b)
{
	 var tr, tg, tb;

	for (var v=0; v<this.layerSizeX*this.layerSizeY; v++)
	{
		tr=this.layers[src][v].r+r;
		tg=this.layers[src][v].g+g;
		tb=this.layers[src][v].b+b;

		if (tr>=255) this.layers[dest][v].r=255; else if (tr<=0) this.layers[dest][v].r=0; else this.layers[dest][v].r=tr;
		if (tg>=255) this.layers[dest][v].g=255; else if (tg<=0) this.layers[dest][v].g=0; else this.layers[dest][v].g=tg;
		if (tb>=255) this.layers[dest][v].b=255; else if (tb<=0) this.layers[dest][v].b=0; else this.layers[dest][v].b=tb;
	}
}


Chromanin.prototype.adjustLayerHSV= function adjustLayerHSV( src,  dest,  h,  s,  v)
{
	var th, ts, tv, tr, tg, tb;

	for ( k=0; k<this.layerSizeX*this.layerSizeY; k++)
	{
		rgb2hsv(this.layers[src][k].r, this.layers[src][k].g, this.layers[src][k].b, th, ts, tv);
		th+=h;
		ts+=s;
		tv+=v;
		if (ts>1) ts=1; else if (ts<0) ts=0;
		if (tv>255) tv=255; else if (tv<0) tv=0;
		hsv2rgb(th, ts, tv, tr, tg, tb);
		this.layers[dest][k].r=getlo(tr);
		this.layers[dest][k].g=getlo(tg);
		this.layers[dest][k].b=getlo(tb);
	}
}


Chromanin.prototype.sineLayerRGB= function sineLayerRGB( src,  dest,  r,  g,  b)
{
	r*=Math.PI;
	g*=Math.PI;
	b*=Math.PI;
	for (var v=0; v<this.layerSizeX*this.layerSizeY; v++)
	{
		this.layers[dest][v].r=(127.5*(Math.sin(r*this.layers[src][v].r)+1));
		this.layers[dest][v].g=(127.5*(Math.sin(g*this.layers[src][v].g)+1));
		this.layers[dest][v].b=(127.5*(Math.sin(b*this.layers[src][v].b)+1));
	}
}


Chromanin.prototype.equalizeRGB= function equalizeRGB( src,  dest)
{
	var histogramR=new Array(256);
	var histogramG=new Array(256);
	var histogramB=new Array(256);
	var sumR, sumG, sumB, pDiv;

	//memset(histogramR, 0, sizeof()*256);
	//memset(histogramG, 0, sizeof()*256);
	//memset(histogramB, 0, sizeof()*256);

	for (var v=0; v<this.layerSizeX*this.layerSizeY; v++)
	{
		histogramR[this.layers[src][v].r]++;
		histogramG[this.layers[src][v].g]++;
		histogramB[this.layers[src][v].b]++;
	}

	sumR=sumG=sumB=0;
	pDiv=255/(this.layerSizeX*this.layerSizeY);
	for (v=0; v<256; v++)			// v is already defined by for (var v=0; v<this.layerSizeX*this.layerSizeY; v++)!!!
	{
		sumR+=histogramR[v]*pDiv;
		histogramR[v]=sumR;
		sumG+=histogramG[v]*pDiv;
		histogramG[v]=sumG;
		sumB+=histogramB[v]*pDiv;
		histogramB[v]=sumB;
	}

	for (v=0; v<this.layerSizeX*this.layerSizeY; v++)
	{
		this.layers[dest][v].r=histogramR[this.layers[src][v].r];
		this.layers[dest][v].g=histogramG[this.layers[src][v].g];
		this.layers[dest][v].b=histogramB[this.layers[src][v].b];
	}
}


Chromanin.prototype.stretchRGB= function stretchRGB( src,  dest)
{
	var histogramR=new Array(256);
	var histogramG=new Array(256);
	var histogramB=new Array(256);
	var sumR, sumG, sumB, pDivR, pDivG, pDivB;
	var minR, minG, minB, maxR, maxG, maxB;

	//memset(histogramR, 0, sizeof()*256);
	//memset(histogramG, 0, sizeof()*256);
	//memset(histogramB, 0, sizeof()*256);

	for (var v=0; v<this.layerSizeX*this.layerSizeY; v++)
	{
		histogramR[this.layers[src][v].r]++;
		histogramG[this.layers[src][v].g]++;
		histogramB[this.layers[src][v].b]++;
	}

	minR=minG=minB=-1;
	maxR=maxG=maxB=0;
	for (v=0; v<256; v++)
	{
		if (histogramR[v]!=0){maxR=v; if(minR==-1)minR=v;}
		if (histogramG[v]!=0){maxG=v; if(minG==-1)minG=v;}
		if (histogramB[v]!=0){maxB=v; if(minB==-1)minB=v;}
	}

	sumR=minR; sumG=minG; sumB=minB;
	pDivR=(maxR-minR)/(this.layerSizeX*this.layerSizeY);
	pDivG=(maxG-minG)/(this.layerSizeX*this.layerSizeY);
	pDivB=(maxB-minB)/(this.layerSizeX*this.layerSizeY);
	for (v=0; v<256; v++)			// v is already defined by for (var v=0; v<this.layerSizeX*this.layerSizeY; v++)!!!
	{
		sumR+=histogramR[v]*pDivR;
		histogramR[v]=sumR;
		sumG+=histogramG[v]*pDivG;
		histogramG[v]=sumG;
		sumB+=histogramB[v]*pDivB;
		histogramB[v]=sumB;
	}

	for (v=0; v<this.layerSizeX*this.layerSizeY; v++)
	{
		this.layers[dest][v].r=histogramR[this.layers[src][v].r];
		this.layers[dest][v].g=histogramG[this.layers[src][v].g];
		this.layers[dest][v].b=histogramB[this.layers[src][v].b];
	}
}


Chromanin.prototype.invertLayer= function invertLayer( src,  dest)
{
	for (var v=0; v<this.layerSizeX*this.layerSizeY; v++)
	{
		this.layers[dest][v].r=255-this.layers[src][v].r;//~this.layers[src][v].r;
		this.layers[dest][v].g=255-this.layers[src][v].g;//~this.layers[src][v].g;
		this.layers[dest][v].b=255-this.layers[src][v].b;//~this.layers[src][v].b;
		//this.layers[dest][v].a=~this.layers[src][v].a;
		}
}


Chromanin.prototype.matrixLayer= function matrixLayer( src,  dest,  myabs,  matrix)
{
	var r, g, b;
	var offset;

	this.copyTemp(src);

	for (var y=0; y<this.layerSizeY; y++)
		for (var x=0; x<this.layerSizeX; x++)
		{
			r=g=b=0;
			for (var my=0; my<3; my++)
				for (var mx=0; mx<3; mx++)
				{
					offset=((x-1+mx)&this.andLayerSizeX)+((y-1+my)&this.andLayerSizeY)*this.layerSizeX;
					r+=matrix[mx][my]*this.layers[this.TEMPL][offset].r;
					g+=matrix[mx][my]*this.layers[this.TEMPL][offset].g;
					b+=matrix[mx][my]*this.layers[this.TEMPL][offset].b;
				}
			if (myabs==true)
			{
				r=Math.abs(r);
				g=Math.abs(g);
				b=Math.abs(b);
			}

			if (r<0) r=0; else if (r>255) r=255;
			if (g<0) g=0; else if (g>255) g=255;
			if (b<0) b=0; else if (b>255) b=255;

			this.layers[dest][x+y*this.layerSizeX].r=r;
			this.layers[dest][x+y*this.layerSizeX].g=g;
			this.layers[dest][x+y*this.layerSizeX].b=b;
		}
}


Chromanin.prototype.embossLayer= function embossLayer( src,  dest)
{
	var r1, g1, b1, r2, g2, b2;
	var offset, offsetxm1, offsetxp1, offsetym1, offsetyp1;

	this.copyTemp(src);

	for (var y=0; y<this.layerSizeY; y++)
	{
		offsetym1=(getlo(y-1)&this.andLayerSizeY)*this.layerSizeX;
		offset=y*this.layerSizeX;
		offsetyp1=(getlo(y+1)&this.andLayerSizeY)*this.layerSizeX;
		for (var x=0; x<this.layerSizeX; x++)
		{
			offsetxm1=(getlo(x-1)&this.andLayerSizeX);
			offsetxp1=(getlo(x+1)&this.andLayerSizeX);
			r1=128
				-this.layers[this.TEMPL][offsetxm1+offsetym1].r
				-this.layers[this.TEMPL][offsetxm1+offset].r
				-this.layers[this.TEMPL][offsetxm1+offsetyp1].r
				+this.layers[this.TEMPL][offsetxp1+offsetym1].r
				+this.layers[this.TEMPL][offsetxp1+offset].r
				+this.layers[this.TEMPL][offsetxp1+offsetyp1].r;
			g1=128
				-this.layers[this.TEMPL][offsetxm1+offsetym1].g
				-this.layers[this.TEMPL][offsetxm1+offset].g
				-this.layers[this.TEMPL][offsetxm1+offsetyp1].g
				+this.layers[this.TEMPL][offsetxp1+offsetym1].g
				+this.layers[this.TEMPL][offsetxp1+offset].g
				+this.layers[this.TEMPL][offsetxp1+offsetyp1].g;
			b1=128
				-this.layers[this.TEMPL][offsetxm1+offsetym1].b
				-this.layers[this.TEMPL][offsetxm1+offset].b
				-this.layers[this.TEMPL][offsetxm1+offsetyp1].b
				+this.layers[this.TEMPL][offsetxp1+offsetym1].b
				+this.layers[this.TEMPL][offsetxp1+offset].b
				+this.layers[this.TEMPL][offsetxp1+offsetyp1].b;
			r2=128
				-this.layers[this.TEMPL][offsetym1+offsetxm1].r
				-this.layers[this.TEMPL][offsetym1+x].r
				-this.layers[this.TEMPL][offsetym1+offsetxp1].r
				+this.layers[this.TEMPL][offsetyp1+offsetxm1].r
				+this.layers[this.TEMPL][offsetyp1+x].r
				+this.layers[this.TEMPL][offsetyp1+offsetxp1].r;
			g2=128
				-this.layers[this.TEMPL][offsetym1+offsetxm1].g
				-this.layers[this.TEMPL][offsetym1+x].g
				-this.layers[this.TEMPL][offsetym1+offsetxp1].g
				+this.layers[this.TEMPL][offsetyp1+offsetxm1].g
				+this.layers[this.TEMPL][offsetyp1+x].g
				+this.layers[this.TEMPL][offsetyp1+offsetxp1].g;
			b2=128
				-this.layers[this.TEMPL][offsetym1+offsetxm1].b
				-this.layers[this.TEMPL][offsetym1+x].b
				-this.layers[this.TEMPL][offsetym1+offsetxp1].b
				+this.layers[this.TEMPL][offsetyp1+offsetxm1].b
				+this.layers[this.TEMPL][offsetyp1+x].b
				+this.layers[this.TEMPL][offsetyp1+offsetxp1].b;
			r1=Math.sqrt((r1*r1+r2*r2));
			g1=Math.sqrt((g1*g1+g2*g2));
			b1=Math.sqrt((b1*b1+b2*b2));
			if (r1>255) r1=255;
			if (g1>255) g1=255;
			if (b1>255) b1=255;
			this.layers[dest][x+offset].r=r1;
			this.layers[dest][x+offset].g=g1;
			this.layers[dest][x+offset].b=b1;
		}
	}
}


// table used for shift operations
powtab = [ 1, 2, 4, 8, 16, 32, 64, 128, 256, 512];

function mylshift(a,shift){
// logical shift left
if (shift > 7) {
     a = 0; // if shifting more than 15 bits to the left, value is always zero
} else {
     a *= powtab[shift];
}
return a;
}

function myrshift(a,shift){
// logical shift right (unsigned)
if (shift > 7) {
    a = 0; // more than 15, becomes zero
} else {//if (shift > 0) {
    //if (a < 0) {
        // deal with the sign bit (15)
      //  a += -32768;
       // a /= powtab[shift];
       // a += powtab[15 - shift];
    //} else {
        a /= powtab[shift];
    //}
}
return a;
}



function myor(a,b)
{
// OR (|)
c = 0;
for (x = 0; x <= 7; ++x) {
    c += c;
    if (a < 0) {
        c += 1;
    } else if (b < 0) {
        c += 1;
    }
    a += a;
    b += b;
}
return c;
}

function gethi(i)
{
    lo = i % 256;
    hi = (i-lo)/256;
	return hi;
}

function getlo(i)
{
return i % 256;
}
Chromanin.prototype.woodLayer= function woodLayer( src,  dest,  b)
{
	var bm8=8-b;
	for (var v=0; v<this.layerSizeX*this.layerSizeY; v++)
	{

		//this.layers[dest][v].r=(this.layers[src][v].r<<b)|(this.layers[src][v].r>>bm8);
		//this.layers[dest][v].g=(this.layers[src][v].g<<b)|(this.layers[src][v].g>>bm8);
		//this.layers[dest][v].b=(this.layers[src][v].b<<b)|(this.layers[src][v].b>>bm8);
		
		this.layers[dest][v].r=getlo((mylshift(this.layers[src][v].r,b))|(myrshift(this.layers[src][v].r,bm8)));
		this.layers[dest][v].g=getlo((mylshift(this.layers[src][v].g,b))|(myrshift(this.layers[src][v].g,bm8)));
		this.layers[dest][v].b=getlo((mylshift(this.layers[src][v].b,b))|(myrshift(this.layers[src][v].b,bm8)));
		
	}
}


Chromanin.prototype.blurLayer= function blurLayer( src,  dest)
{

	var matrix=new Array(3);
	matrix[0]=new Array(3);
	matrix[1]=new Array(3);
	matrix[2]=new Array(3);
	
	matrix[0][0]=0.0625; matrix[1][0]=0.125; matrix[2][0]=0.0625;
	matrix[0][1]=0.125; matrix[1][1]=0.25; matrix[2][1]=0.125;
	matrix[0][2]=0.0625; matrix[1][2]=0.125; matrix[2][2]=0.0625;
	this.matrixLayer(src, dest, false, matrix);
}


Chromanin.prototype.edgeHLayer= function edgeHLayer( src,  dest)
{
	var matrix=new Array(3);
	matrix[0]=new Array(3);
	matrix[1]=new Array(3);
	matrix[2]=new Array(3);
	
	matrix[0][0]=2; matrix[1][0]=4; matrix[2][0]=2;
	matrix[0][1]=0; matrix[1][1]=0; matrix[2][1]=0;
	matrix[0][2]=-2; matrix[1][2]=-4; matrix[2][2]=-2;
	this.matrixLayer(src, dest, true, matrix);
}


Chromanin.prototype.edgeVLayer= function edgeVLayer( src,  dest)
{
	var matrix=new Array(3);
	matrix[0]=new Array(3);
	matrix[1]=new Array(3);
	matrix[2]=new Array(3);
	
	matrix[0][0]=2; matrix[1][0]=0; matrix[2][0]=-2;
	matrix[0][1]=4; matrix[1][1]=0; matrix[2][1]=-4;
	matrix[0][2]=2; matrix[1][2]=0; matrix[2][2]=-2;
	this.matrixLayer(src, dest, true, matrix);
}


Chromanin.prototype.sharpenLayer= function sharpenLayer( src,  dest)
{
	var matrix=new Array(3);
	matrix[0]=new Array(3);
	matrix[1]=new Array(3);
	matrix[2]=new Array(3);
	
	matrix[0][0]=-0.125; matrix[1][0]=-0.25; matrix[2][0]=-0.125;
	matrix[0][1]=-0.25; matrix[1][1]=2.5; matrix[2][1]=-0.25;
	matrix[0][2]=-0.125; matrix[1][2]=-0.25; matrix[2][2]=-0.125;
	this.matrixLayer(src, dest, false, matrix);
}


Chromanin.prototype.motionBlur= function motionBlur( src,  dest,  s)
{
	var sq, ts, r, g, b, offset, offset2;

	this.copyTemp(src);
	sq=(s+1)*(s+1);
	for (var y=0; y<this.layerSizeY; y++)
	{
		offset=y*this.layerSizeX;
		for (var x=0; x<this.layerSizeX; x++)
		{
			r=g=b=0;
			for ( t=-s; t<=s; t++)
			{
				offset2=offset+(getlo(x+t)&this.andLayerSizeX);
				ts=(s+1-Math.abs(t));
				r+=this.layers[this.TEMPL][offset2].r*ts;
				g+=this.layers[this.TEMPL][offset2].g*ts;
				b+=this.layers[this.TEMPL][offset2].b*ts;
			}
			this.layers[dest][offset+x].r=(r/sq);
			this.layers[dest][offset+x].g=(g/sq);
			this.layers[dest][offset+x].b=(b/sq);
		}
	}
}


Chromanin.prototype.makeTilable= function makeTilable( src,  dest,  s)
{
	var offset, offset2, sx, sy;
	var sq, sr, sd, srm;

	this.copyTemp(src);
	s=this.layerSizeX/2-s;
	sq=(s*s);
	sd=(this.layerSizeX/2)*(this.layerSizeY/2)-sq;
	if (sd!=0)
	   {
		sd=0.75/sd;
		}
	else
	{
		sd=75000;
		}
	for (var y=0; y<this.layerSizeY; y++)
	{
		offset=y*this.layerSizeX;
		offset2=(this.layerSizeY-1-y)*this.layerSizeX;
		sy=y-this.layerSizeY/2;
		sy*=sy;
		for (var x=0; x<this.layerSizeX; x++)
		{
			sx=x-this.layerSizeX/2;
			sr=sx*sx+sy-sq;
			if (sr>0)
			{
				sr*=sd;
				if (sr>0.75)
				{
					sr=0.25;
					srm=0.25;
				} else {
					srm=1-sr;
					sr/=3;
				}
				this.layers[dest][offset+x].r=(this.layers[this.TEMPL][offset+x].r*srm+(this.layers[this.TEMPL][offset+this.layerSizeX-1-x].r+this.layers[this.TEMPL][offset2+this.layerSizeX-1-x].r+this.layers[this.TEMPL][offset2+x].r)*sr);
				this.layers[dest][offset+x].g=(this.layers[this.TEMPL][offset+x].g*srm+(this.layers[this.TEMPL][offset+this.layerSizeX-1-x].g+this.layers[this.TEMPL][offset2+this.layerSizeX-1-x].g+this.layers[this.TEMPL][offset2+x].g)*sr);
				this.layers[dest][offset+x].b=(this.layers[this.TEMPL][offset+x].b*srm+(this.layers[this.TEMPL][offset+this.layerSizeX-1-x].b+this.layers[this.TEMPL][offset2+this.layerSizeX-1-x].b+this.layers[this.TEMPL][offset2+x].b)*sr);
			}
		}
	}
}


function median(v)
{
	 var a= Array(5);
	 var i, j, last;
	//memset(a, 0, 5);
	for (i=0; i<9; i++)
	{
		last=4;
		for (j=4; j>=0; j--)
			if (a[j]<=v[i])
				last=j;
		for (j=3; j>=last; j--)
			a[j+1]=a[j];
		a[last]=v[i];
	}
	return a[4];
}


Chromanin.prototype.medianLayer= function medianLayer( src,  dest)
{
	var offset = new Array(9);
	var colors=new Array(9);
	var i;

	this.copyTemp(src);

	for (var y=0; y<this.layerSizeY; y++)
		for (var x=0; x<this.layerSizeX; x++)
		{
			for (i=0; i<9; i++)
				offset[i]=((x-1+(i%3))&this.andLayerSizeX)+((y-1+(i/3))&this.andLayerSizeY)*this.layerSizeX;
			for (i=0; i<9; i++)
				colors[i]=this.layers[this.TEMPL][offset[i]].r;
			this.layers[dest][offset[4]].r=median(colors);
			for (i=0; i<9; i++)
				colors[i]=this.layers[this.TEMPL][offset[i]].g;
			this.layers[dest][offset[4]].g=median(colors);
			for (i=0; i<9; i++)
				colors[i]=this.layers[this.TEMPL][offset[i]].b;
			this.layers[dest][offset[4]].b=median(colors);
		}
}

Chromanin.prototype.copyTemp= function copyTemp(src)
{

for (var x=0; x<(this.layerSizeX * this.layerSizeY); x+=1)
		{
		    
		    this.layers[5][x].r=this.layers[src][x].r;
			this.layers[5][x].g=this.layers[src][x].g;
			this.layers[5][x].b=this.layers[src][x].b;
			this.layers[5][x].a=this.layers[src][x].a;
			
		}

}

Chromanin.prototype.erodeLayer= function erodeLayer( src,  dest)
{
	var offset, offsetym1, offsetyp1, offsetxm1, offsetxp1;
	var r, g, b;

	
	this.copyTemp(src);
	
	
	for (var y=0; y<this.layerSizeY; y++)
	{
		offset=y*this.layerSizeX;
		offsetym1=getlo((y-1)&this.andLayerSizeY)*this.layerSizeX;
		offsetyp1=getlo((y+1)&this.andLayerSizeY)*this.layerSizeX;
		for (var x=0; x<this.layerSizeX; x++)
		{
			offsetxm1=getlo(x-1)&this.andLayerSizeX;
			offsetxp1=getlo(x+1)&this.andLayerSizeX;
			

			r=	Math.min(this.layers[this.TEMPL][offsetym1+offsetxm1].r,
				Math.min(this.layers[this.TEMPL][offsetym1+x].r,
				Math.min(this.layers[this.TEMPL][offsetym1+offsetxp1].r,
				Math.min(this.layers[this.TEMPL][offset+offsetxm1].r,
				Math.min(this.layers[this.TEMPL][offset+x].r,
				Math.min(this.layers[this.TEMPL][offset+offsetxp1].r,
				Math.min(this.layers[this.TEMPL][offsetyp1+offsetxm1].r,
				Math.min(this.layers[this.TEMPL][offsetyp1+x].r,
				this.layers[this.TEMPL][offsetyp1+offsetxp1].r))))))));
			g=	Math.min(this.layers[this.TEMPL][offsetym1+offsetxm1].g,
				Math.min(this.layers[this.TEMPL][offsetym1+x].g,
				Math.min(this.layers[this.TEMPL][offsetym1+offsetxp1].g,
				Math.min(this.layers[this.TEMPL][offset+offsetxm1].g,
				Math.min(this.layers[this.TEMPL][offset+x].g,
				Math.min(this.layers[this.TEMPL][offset+offsetxp1].g,
				Math.min(this.layers[this.TEMPL][offsetyp1+offsetxm1].g,
				Math.min(this.layers[this.TEMPL][offsetyp1+x].g,
				this.layers[this.TEMPL][offsetyp1+offsetxp1].g))))))));
			b=	Math.min(this.layers[this.TEMPL][offsetym1+offsetxm1].b,
				Math.min(this.layers[this.TEMPL][offsetym1+x].b,
				Math.min(this.layers[this.TEMPL][offsetym1+offsetxp1].b,
				Math.min(this.layers[this.TEMPL][offset+offsetxm1].b,
				Math.min(this.layers[this.TEMPL][offset+x].b,
				Math.min(this.layers[this.TEMPL][offset+offsetxp1].b,
				Math.min(this.layers[this.TEMPL][offsetyp1+offsetxm1].b,
				Math.min(this.layers[this.TEMPL][offsetyp1+x].b,
				this.layers[this.TEMPL][offsetyp1+offsetxp1].b))))))));
				
			this.layers[dest][offset+x].r=r;
			this.layers[dest][offset+x].g=g;
			this.layers[dest][offset+x].b=b;
		}
	}
}


Chromanin.prototype.dilateLayer= function dilateLayer( src,  dest)
{
	var offset, offsetym1, offsetyp1, offsetxm1, offsetxp1;
	var r, g, b;

	this.copyTemp(src);
	
	for (var y=0; y<this.layerSizeY; y++)
	{
		offset=y*this.layerSizeX;
		offsetym1=(getlo(y-1)&this.andLayerSizeY)*this.layerSizeX;
		offsetyp1=(getlo(y+1)&this.andLayerSizeY)*this.layerSizeX;
		for (var x=0; x<this.layerSizeX; x++)
		{
			offsetxm1=getlo(x-1)&this.andLayerSizeX;
			offsetxp1=getlo(x+1)&this.andLayerSizeX;

			r=	Math.max(this.layers[this.TEMPL][offsetym1+offsetxm1].r,
				Math.max(this.layers[this.TEMPL][offsetym1+x].r,
				Math.max(this.layers[this.TEMPL][offsetym1+offsetxp1].r,
				Math.max(this.layers[this.TEMPL][offset+offsetxm1].r,
				Math.max(this.layers[this.TEMPL][offset+x].r,
				Math.max(this.layers[this.TEMPL][offset+offsetxp1].r,
				Math.max(this.layers[this.TEMPL][offsetyp1+offsetxm1].r,
				Math.max(this.layers[this.TEMPL][offsetyp1+x].r,
				this.layers[this.TEMPL][offsetyp1+offsetxp1].r))))))));
			g=	Math.max(this.layers[this.TEMPL][offsetym1+offsetxm1].g,
				Math.max(this.layers[this.TEMPL][offsetym1+x].g,
				Math.max(this.layers[this.TEMPL][offsetym1+offsetxp1].g,
				Math.max(this.layers[this.TEMPL][offset+offsetxm1].g,
				Math.max(this.layers[this.TEMPL][offset+x].g,
				Math.max(this.layers[this.TEMPL][offset+offsetxp1].g,
				Math.max(this.layers[this.TEMPL][offsetyp1+offsetxm1].g,
				Math.max(this.layers[this.TEMPL][offsetyp1+x].g,
				this.layers[this.TEMPL][offsetyp1+offsetxp1].g))))))));
			b=	Math.max(this.layers[this.TEMPL][offsetym1+offsetxm1].b,
				Math.max(this.layers[this.TEMPL][offsetym1+x].b,
				Math.max(this.layers[this.TEMPL][offsetym1+offsetxp1].b,
				Math.max(this.layers[this.TEMPL][offset+offsetxm1].b,
				Math.max(this.layers[this.TEMPL][offset+x].b,
				Math.max(this.layers[this.TEMPL][offset+offsetxp1].b,
				Math.max(this.layers[this.TEMPL][offsetyp1+offsetxm1].b,
				Math.max(this.layers[this.TEMPL][offsetyp1+x].b,
				this.layers[this.TEMPL][offsetyp1+offsetxp1].b))))))));
			this.layers[dest][offset+x].r=r;
			this.layers[dest][offset+x].g=g;
			this.layers[dest][offset+x].b=b;
		}
	}
}


Chromanin.prototype.sineDistort= function sineDistort( src,  dest,  dx,  depthX,  dy,  depthY)
{
	this.copyTemp(src);
	for (var y=0; y<this.layerSizeY; y++)
		for (var x=0; x<this.layerSizeX; x++)
			this.layers[dest][x+y*this.layerSizeX]=this.getBilerpixel(this.TEMPL, Math.sin(dx*y)*depthX+x, Math.sin(dy*x)*depthY+y);
}


Chromanin.prototype.twirlLayer= function twirlLayer( src,  dest,  rot,  scale)
{
	 var ooscale, a, b, d, winkel, cw, sw, na, nb;
	 var ina, inb, inbp, inap1, inbp1;
	var corner = [];
	 corner[0] =new RGBA();
	 corner[1] =new RGBA();
	 corner[2] =new RGBA();
	 corner[3] =new RGBA();

	this.copyTemp(src);

	ooscale=1/(scale*Math.sqrt(2*this.layerSizeX*this.layerSizeY));
	for (var y=0; y<this.layerSizeY; y++)
	{
		b=(y-this.layerSizeY/2);
		for (var x=0; x<this.layerSizeX; x++)
		{
			a=(x-this.layerSizeX/2);
			d=-Math.sqrt(a*a+b*b)+this.layerSizeX/2;
			if (d<=0)
			{
				na=(ina=x);
				nb=(inb=y);
			} else {
				winkel=rot*d*d*ooscale;
				cw=Math.cos(winkel);
				sw=Math.sin(winkel);
				na=a*cw-b*sw+this.layerSizeX/2;
				nb=a*sw+b*cw+this.layerSizeY/2;
				ina=na;
				inb=nb;
			}
			inbp=(inb&this.andLayerSizeY)*this.layerSizeX;
			inbp1=((inb+1)&this.andLayerSizeY)*this.layerSizeX;
			inap1=(ina+1)&this.andLayerSizeX;
			corner[0]=this.layers[this.TEMPL][(ina&this.andLayerSizeX)+inbp];
			corner[1]=this.layers[this.TEMPL][inap1+inbp];
			corner[2]=this.layers[this.TEMPL][(ina&this.andLayerSizeX)+inbp1];
			corner[3]=this.layers[this.TEMPL][inap1+inbp1];
			this.layers[dest][x+y*this.layerSizeX]=cosineInterpolate(corner, na-ina, nb-inb);
		}
	}
}


Chromanin.prototype.tileLayer= function tileLayer( src,  dest)
{
	 var offset, offset2, offset3;

	this.copyTemp(src);

	for (var y=0; y<this.layerSizeY; y++)
	{
		offset=y*this.layerSizeX;
		offset2=((y*2)&this.andLayerSizeY)*this.layerSizeX;
		for (var x=0; x<this.layerSizeX; x++)
		{
			offset3=((x*2)&this.andLayerSizeX)+offset2;
			this.layers[dest][x+offset].r=(this.layers[this.TEMPL][offset3].r+this.layers[this.TEMPL][offset3+1].r+this.layers[this.TEMPL][offset3+this.layerSizeX].r+this.layers[this.TEMPL][offset3+this.layerSizeX+1].r)>>2;
			this.layers[dest][x+offset].g=(this.layers[this.TEMPL][offset3].g+this.layers[this.TEMPL][offset3+1].g+this.layers[this.TEMPL][offset3+this.layerSizeX].g+this.layers[this.TEMPL][offset3+this.layerSizeX+1].g)>>2;
			this.layers[dest][x+offset].b=(this.layers[this.TEMPL][offset3].b+this.layers[this.TEMPL][offset3+1].b+this.layers[this.TEMPL][offset3+this.layerSizeX].b+this.layers[this.TEMPL][offset3+this.layerSizeX+1].b)>>2;
		}
	}
}


Chromanin.prototype.noiseDistort= function noiseDistort( src,  dest,  seed,  radius)
{
	var dx, dy;
	this.copyTemp(src);
	this.seedValue=seed;
	radius=16-radius;
	for (var y=0; y<this.layerSizeY; y++)
		for (var x=0; x<this.layerSizeX; x++)
		{
			
			dx=myRandom()+radius;//myrshift((myRandom()),radius);
			dy=myRandom()+radius;//myrshift((myRandom()),radius);
			this.layers[dest][x+y*this.layerSizeX].r=this.layers[this.TEMPL][(getlo(x+dx)&this.andLayerSizeX)+(getlo(y+dy)&this.andLayerSizeY)*this.layerSizeX].r;
			this.layers[dest][x+y*this.layerSizeX].g=this.layers[this.TEMPL][(getlo(x+dx)&this.andLayerSizeX)+(getlo(y+dy)&this.andLayerSizeY)*this.layerSizeX].g;
			this.layers[dest][x+y*this.layerSizeX].b=this.layers[this.TEMPL][(getlo(x+dx)&this.andLayerSizeX)+(getlo(y+dy)&this.andLayerSizeY)*this.layerSizeX].b;
	
		}
}


Chromanin.prototype.moveDistort= function moveDistort( src,  dest,  dx,  dy)
{
	this.copyTemp(src);
	for (var y=0; y<this.layerSizeY; y++)
		for (var x=0; x<this.layerSizeX; x++)
		{
			this.layers[dest][x+y*this.layerSizeX].r=this.layers[this.TEMPL][(getlo(x+dx)&this.andLayerSizeX)+(getlo(y+dy)&this.andLayerSizeY)*this.layerSizeX].r;
			this.layers[dest][x+y*this.layerSizeX].g=this.layers[this.TEMPL][(getlo(x+dx)&this.andLayerSizeX)+(getlo(y+dy)&this.andLayerSizeY)*this.layerSizeX].g;
			this.layers[dest][x+y*this.layerSizeX].b=this.layers[this.TEMPL][(getlo(x+dx)&this.andLayerSizeX)+(getlo(y+dy)&this.andLayerSizeY)*this.layerSizeX].b;
			}
}


Chromanin.prototype.move= function move(src,inds,dest,indd,size)
{

for (var x=0; x<(size)-1; x+=1)
		{
		    //alert(src);
			//alert(dest);
			//alert(inds);
			//alert(indd);
		    this.layers[dest][indd].r=this.layers[src][inds].r;
			this.layers[dest][indd].g=this.layers[src][inds].g;
			this.layers[dest][indd].b=this.layers[src][inds].b;
			this.layers[dest][indd].a=this.layers[src][inds].a;
			indd++;
			inds++;
		}

}


Chromanin.prototype.MirrorCorner= function MirrorCorner(c0,dest)
 {
  var xc,yc,offset;
 
   switch(c0){

   case 0: 
         for (var yc=0;yc<(this.layerSizeY/2)-1;yc++)
		   {
          offset=yc*this.layerSizeX;
          for (var xc=0;xc<(this.layerSizeX/2)-1;xc++)
		  {
             this.layers[dest][Math.floor(offset+(this.layerSizeX-xc))].r=this.layers[dest][offset+xc].r;
			 this.layers[dest][Math.floor(offset+(this.layerSizeX-xc))].g=this.layers[dest][offset+xc].g;
			 this.layers[dest][Math.floor(offset+(this.layerSizeX-xc))].b=this.layers[dest][offset+xc].b;
		   }
		 
		 }
       break;
    case 1: 
	     for (var yc=0;yc<(this.layerSizeY/2)-1;yc++)
		 {
         this.move(dest,yc*this.layerSizeX+(this.layerSizeX/2),dest,(this.layerSizeY-yc-1)*this.layerSizeX+(this.layerSizeX/2),(this.layerSizeX/2));
		  }
		break;
	case 2: 
	   for (var yc=0;yc<(this.layerSizeY/2)-1;yc++)
		 {
        offset=(yc+(this.layerSizeY/2))*this.layerSizeX;
        for (var xc=0;xc<(this.layerSizeX/2)-1;xc++)
		 {
		 //alert(offset+xc);
         this.layers[dest][Math.floor(offset+xc)].r=this.layers[dest][offset+(this.layerSizeX-xc)].r;
		 this.layers[dest][Math.floor(offset+xc)].g=this.layers[dest][offset+(this.layerSizeX-xc)].g;
		 this.layers[dest][Math.floor(offset+xc)].b=this.layers[dest][offset+(this.layerSizeX-xc)].b;
		 }
		}
       break;
    case 3: 

		for (var yc=0;yc<(this.layerSizeY/2)-1;yc++)
		{
		 this.move(dest,(this.layerSizeY-yc-1)*this.layerSizeX,dest,yc*this.layerSizeX,(this.layerSizeX/2));
		 }
		break;
  }
 }


Chromanin.prototype.kaleidLayer= function kaleidLayer( src,  dest,  corner)
{
	
 var y;

 

  corner=corner-1;

  for (var y=0; y<(this.layerSizeY/2);y++)
  {
   //this.move(src,(y+(corner/2)*(this.layerSizeY/2))*this.layerSizeX+(corner%2)*(this.layerSizeX/2),dest,(y+(corner/2)*(this.layerSizeY/2))*this.layerSizeX+(corner%2)*(this.layerSizeX/2),(this.layerSizeX/2));

  this.MirrorCorner(corner,dest);
  this.MirrorCorner((corner+1)%4,dest);
  this.MirrorCorner((corner+2)%4,dest);
}
	
}


Chromanin.prototype.tunnelDistort= function tunnelDistort( src,  dest,  f)
{
	 var ina, inb, inap1, inbp, inbp1;
	 var a, b, na, nb;
	 var corner = [];
	 corner[0] =new RGBA();
	 corner[1] =new RGBA();
	 corner[2] =new RGBA();
	 corner[3] =new RGBA();
	 
	 
	 var arct, lsd2p;

	this.copyTemp(src);
	var lsd2p=this.layerSizeX/(2*Math.PI);
	for (var y=0; y<this.layerSizeY; y++)
	{
		b=-0.5*this.layerSizeY+y;
		for (var x=0; x<this.layerSizeX; x++)
		{
			a=-0.5*this.layerSizeX+x;
			if (a!=0)
			{
				a=1/a;
				arct=Math.atan(b*a);
				if (a>0)
					na=lsd2p*arct+this.layerSizeX/2;
				else
					na=lsd2p*arct;
				nb=Math.abs(Math.cos(arct)*f*a);
				ina=na;
				inb=nb;
				inap1=(ina+1)&this.andLayerSizeX;
				inbp=(inb&this.andLayerSizeY)*this.layerSizeX;
				inbp1=((inb+1)&this.andLayerSizeY)*this.layerSizeX;
				corner[0]=this.layers[this.TEMPL][(ina&this.andLayerSizeX)+inbp];
				corner[1]=this.layers[this.TEMPL][inap1+inbp];
				corner[2]=this.layers[this.TEMPL][(ina&this.andLayerSizeX)+inbp1];
				corner[3]=this.layers[this.TEMPL][inap1+inbp1];
			}
			this.layers[dest][x+y*this.layerSizeX]=cosineInterpolate(corner, na-ina, nb-inb);
		} 
	}
}


Chromanin.prototype.sculptureLayer= function sculptureLayer( src,  dest)
{
	var ipi=255.0/(2.0*3.1415926536);
    var x,y,offset,offsetym1,offsetyp1,offsetxm1,offsetxp1;
    var r1,r2,g1,g2,b1,b2;
    var a;

	this.copyTemp(src);

  for (var y=0; y<this.layerSizeY; y++)
  {
   offset=y*this.layerSizeX;
   offsetym1=((y-1)&this.andLayerSizeY)*this.layerSizeX;
   offsetyp1=((y+1)&this.andLayerSizeY)*this.layerSizeX;
   for (var x=0;x<this.layerSizeX; x++)
   {
    offsetxm1=((x-1)&this.andLayerSizeX);
    offsetxp1=((x+1)&this.andLayerSizeX);

	
    r1= this.layers[this.TEMPL][offsetxm1 + offsetym1].r
         +2*this.layers[this.TEMPL][offsetxm1 + offset].r
         +this.layers[this.TEMPL][offsetxm1 + offsetyp1].r
         -this.layers[this.TEMPL][offsetxp1 + offsetym1].r
         -2*this.layers[this.TEMPL][offsetxp1 + offset].r
         -this.layers[this.TEMPL][offsetxp1 + offsetyp1].r;
    r2= this.layers[this.TEMPL][offsetym1 + offsetxm1].r
         +2*this.layers[this.TEMPL][offsetym1 + x].r
         +this.layers[this.TEMPL][offsetym1 + offsetxp1].r
         -this.layers[this.TEMPL][offsetyp1 + offsetxm1].r
         -2*this.layers[this.TEMPL][offsetyp1 + x].r
         -this.layers[this.TEMPL][offsetyp1 + offsetxp1].r;

    g1= this.layers[this.TEMPL][offsetxm1 + offsetym1].g
         +2*this.layers[this.TEMPL][offsetxm1 + offset].g
         +this.layers[this.TEMPL][offsetxm1 + offsetyp1].g
         -this.layers[this.TEMPL][offsetxp1 + offsetym1].g
         -2*this.layers[this.TEMPL][offsetxp1 + offset].g
         -this.layers[this.TEMPL][offsetxp1 + offsetyp1].g;
    g2= this.layers[this.TEMPL][offsetym1 + offsetxm1].g
         +2*this.layers[this.TEMPL][offsetym1 + x].g
         +this.layers[this.TEMPL][offsetym1 + offsetxp1].g
         -this.layers[this.TEMPL][offsetyp1 + offsetxm1].g
         -2*this.layers[this.TEMPL][offsetyp1 + x].g
         -this.layers[this.TEMPL][offsetyp1 + offsetxp1].g;

    b1= this.layers[this.TEMPL][offsetxm1 + offsetym1].b
         +2*this.layers[this.TEMPL][offsetxm1 + offset].b
         +this.layers[this.TEMPL][offsetxm1 + offsetyp1].b
         -this.layers[this.TEMPL][offsetxp1 + offsetym1].b
         -2*this.layers[this.TEMPL][offsetxp1 + offset].b
         -this.layers[this.TEMPL][offsetxp1 + offsetyp1].b;
    b2= this.layers[this.TEMPL][offsetym1 + offsetxm1].b
         +2*this.layers[this.TEMPL][offsetym1 + x].b
         +this.layers[this.TEMPL][offsetym1 + offsetxp1].b
         -this.layers[this.TEMPL][offsetyp1 + offsetxm1].b
         -2*this.layers[this.TEMPL][offsetyp1 + x].b
         -this.layers[this.TEMPL][offsetyp1 + offsetxp1].b;

    if (r1==0) {
     if (r2>0) this.layers[dest][x + offset].r=196;
      else if (r2==0) this.layers[dest][x + offset].r=128;
       else this.layers[dest][x + offset].r=64;
    }
    else {
     a=Math.atan(r2/r1);
     if (r1>0) this.layers[dest][x + offset].r=trunc(a*ipi+127.5)
      else this.layers[dest][x + offset].r=trunc(a*ipi);
    }

    if (g1==0) {
     if (g2>0) this.layers[dest][x + offset].g=196;
      else if (g2==0) this.layers[dest][x + offset].g=128;
       else this.layers[dest][x + offset].g=64;
   }
    else {
     a=Math.atan(g2/g1);
     if (g1>0) this.layers[dest][x + offset].g=trunc(a*ipi+127.5);
      else this.layers[dest][x + offset].g=trunc(a*ipi);
    }

    if (b1==0) {
     if (b2>0) this.layers[dest][x + offset].b=196;
      else if (b2==0) this.layers[dest][x + offset].b=128;
       else this.layers[dest][x + offset].b=64;
    }
    else {
     a=Math.atan(b2/b1);
     if (b1>0) this.layers[dest][x + offset].b=trunc(a*ipi+127.5)
      else this.layers[dest][x + offset].b=trunc(a*ipi);
     }
	} 
  }

}

function trunc(n){
   return ~~n;
}
Chromanin.prototype.mapDistort= function mapDistort( src,  dist,  dest,  xd,  yd)
{
	var offset;
	var v;

	this.copyTemp(src);

	for (var y=0; y<this.layerSizeY; y++)
		for (var x=0; x<this.layerSizeX; x++)
		{
			offset=y*this.layerSizeX+x;
			v=Math.max(this.layers[dist][offset].r, Math.max(this.layers[dist][offset].g, this.layers[dist][offset].b));
			this.layers[dest][offset]=this.getBilerpixel(this.TEMPL, xd*v+x, yd*v+y);
		}
}


Chromanin.prototype.addLayers= function addLayers( src1,  src2,  dest,  perc1,  perc2)
{
	 var r, g, b;

	for (var v=0; v<this.layerSizeX*this.layerSizeY; v++)
	{
		r=this.layers[src1][v].r*perc1+this.layers[src2][v].r*perc2;
		g=this.layers[src1][v].g*perc1+this.layers[src2][v].g*perc2;
		b=this.layers[src1][v].b*perc1+this.layers[src2][v].b*perc2;
		if (r>255) r=255; else if (r<0) r=0;
		if (g>255) g=255; else if (g<0) g=0;
		if (b>255) b=255; else if (b<0) b=0;
		this.layers[dest][v].r=r;
		this.layers[dest][v].g=g;
		this.layers[dest][v].b=b;
	}
}


Chromanin.prototype.mulLayers= function mulLayers( src1,  src2,  dest,  perc1,  perc2)
{
	 var r, g, b, perc;

	perc=perc1*perc2/255;
	for (var v=0; v<this.layerSizeX*this.layerSizeY; v++)
	{
		r=this.layers[src1][v].r*this.layers[src2][v].r*perc;
		g=this.layers[src1][v].g*this.layers[src2][v].g*perc;
		b=this.layers[src1][v].b*this.layers[src2][v].b*perc;
		if (r>255) r=255; else if (r<0) r=0;
		if (g>255) g=255; else if (g<0) g=0;
		if (b>255) b=255; else if (b<0) b=0;
		this.layers[dest][v].r=r;
		this.layers[dest][v].g=g;
		this.layers[dest][v].b=b;
	}
}


Chromanin.prototype.xorLayers= function xorLayers( src1,  src2,  dest,  perc1,  perc2)
{
	 var r, g, b;

	for (var v=0; v<this.layerSizeX*this.layerSizeY; v++)
	{
		r=this.layers[src1][v].r^this.layers[src2][v].r;
		g=this.layers[src1][v].g^this.layers[src2][v].g;
		b=this.layers[src1][v].b^this.layers[src2][v].b;
		if (r>255) r=255; else if (r<0) r=0;
		if (g>255) g=255; else if (g<0) g=0;
		if (b>255) b=255; else if (b<0) b=0;
		this.layers[dest][v].r=r;
		this.layers[dest][v].g=g;
		this.layers[dest][v].b=b;
	}
	
}



Chromanin.prototype.andLayers= function andLayers( src1,  src2,  dest,  perc1,  perc2)
{
	var r, g, b;

	for (var v=0; v<this.layerSizeX*this.layerSizeY; v++)
	{
		r=this.layers[src1][v].r&this.layers[src2][v].r;
		g=this.layers[src1][v].g&this.layers[src2][v].g;
		b=this.layers[src1][v].b&this.layers[src2][v].b;
		if (r>255) r=255; else if (r<0) r=0;
		if (g>255) g=255; else if (g<0) g=0;
		if (b>255) b=255; else if (b<0) b=0;
		this.layers[dest][v].r=r;
		this.layers[dest][v].g=g;
		this.layers[dest][v].b=b;
	}
	
}


Chromanin.prototype.orLayers= function orLayers( src1,  src2,  dest,  perc1,  perc2)
{
	var r, g, b;

	for (var v=0; v<this.layerSizeX*this.layerSizeY; v++)
	{
		r=this.layers[src1][v].r|this.layers[src2][v].r;
		g=this.layers[src1][v].g|this.layers[src2][v].g;
		b=this.layers[src1][v].b|this.layers[src2][v].b;
		if (r>255) r=255; else if (r<0) r=0;
		if (g>255) g=255; else if (g<0) g=0;
		if (b>255) b=255; else if (b<0) b=0;
		this.layers[dest][v].r=r;
		this.layers[dest][v].g=g;
		this.layers[dest][v].b=b;
	}
	
}


Chromanin.prototype.randCombineLayers= function randCombineLayers( src1,  src2,  dest,  perc1,  perc2)
{
	for (var v=0; v<this.layerSizeX*this.layerSizeY; v++)
	{
		this.layers[dest][v].r=myRandom()&1?this.layers[src1][v].r:this.layers[src2][v].r;
		this.layers[dest][v].g=myRandom()&1?this.layers[src1][v].g:this.layers[src2][v].g;
		this.layers[dest][v].b=myRandom()&1?this.layers[src1][v].b:this.layers[src2][v].b;
	}
}


Chromanin.prototype.maxCombineLayers= function maxCombineLayers( src1,  src2,  dest,  perc1,  perc2)
{
	for (var v=0; v<this.layerSizeX*this.layerSizeY; v++)
	{
		this.layers[dest][v].r=Math.max(this.layers[src1][v].r, this.layers[src2][v].r);
		this.layers[dest][v].g=Math.max(this.layers[src1][v].g, this.layers[src2][v].g);
		this.layers[dest][v].b=Math.max(this.layers[src1][v].b, this.layers[src2][v].b);
	}
}


Chromanin.prototype.minCombineLayers= function minCombineLayers( src1,  src2,  dest,  perc1,  perc2)
{
	for (var v=0; v<this.layerSizeX*this.layerSizeY; v++)
	{
		this.layers[dest][v].r=Math.min(this.layers[src1][v].r, this.layers[src2][v].r);
		this.layers[dest][v].g=Math.min(this.layers[src1][v].g, this.layers[src2][v].g);
		this.layers[dest][v].b=Math.min(this.layers[src1][v].b, this.layers[src2][v].b);
	}
}

Chromanin.prototype.cellMachine= function cellMachine(l,seed,rule)
{


 var x,y,base_off,m ;
     var c = new RGBA();

  this.seedValue=seed;

  
  c.r=255; c.g=255; c.b=255;

  for (var x=0; x<this.layerSizeX;x++)
  {
   if (((myRandom())>>100)==0) 
   {
     this.layers[l][x]=c;
	 }
}

  base_off=0;
  for (var y=1;y<this.layerSizeY;y++)
  {
   for (var x=0; x<this.layerSizeX;x++)
   {
    if (this.layers[l][((x-1)&this.andLayerSizeX)+base_off].r!=0) 
	  {
	   m=1;
	   }
     else 
	   {
	    m=0;
		}
    if (this.layers[l][x+base_off].r!=0) 
        {	m=m|2;
		 }
		 
    if (this.layers[l][((x+1)&this.andLayerSizeX) +base_off].r!=0) 
	  {
	   m=m|4;
	   }

    if (((1<<m)&rule)!=0)
	 {
	  this.layers[l][x+base_off+this.layerSizeX].r=c.r;
	  this.layers[l][x+base_off+this.layerSizeX].g=c.g;
	  this.layers[l][x+base_off+this.layerSizeX].b=c.b;
	  }
   }
   
   //
   base_off=base_off + this.layerSizeX;
  }

}


