//
//  OPJSupport.cpp
//  OsiriX_Lion
//
//  Created by Aaron Boxer on 1/21/14.
//

#include "OPJSupport.h"
#include "../Binaries/openjpeg/openjpeg.h"

#define J2K_CFMT 0
#define JP2_CFMT 1
#define JPT_CFMT 2


typedef struct decode_info
{
	opj_codec_t *codec;
	opj_stream_t *stream;
	opj_image_t *image;
	OPJ_BOOL   deleteImage;

} decode_info_t;

#define JP2_RFC3745_MAGIC "\x00\x00\x00\x0c\x6a\x50\x20\x20\x0d\x0a\x87\x0a"
#define JP2_MAGIC "\x0d\x0a\x87\x0a"
/* position 45: "\xff\x52" */
#define J2K_CODESTREAM_MAGIC "\xff\x4f\xff\x51"


int buffer_format(opj_buffer_info_t* buf_info)
{
	int magic_format;
	if (!buf_info || buf_info->len < 12)
		return -1;
	if(memcmp(buf_info->buf, JP2_RFC3745_MAGIC, 12) == 0
       || memcmp(buf_info->buf, JP2_MAGIC, 4) == 0)
	{
		magic_format = JP2_CFMT;
	}
	else
	{
		if(memcmp(buf_info->buf, J2K_CODESTREAM_MAGIC, 4) == 0)
		{
			magic_format = J2K_CFMT;
		}
		else
			return -1;
	}
	return magic_format;
}/*  buffer_format() */

const char *clr_space(OPJ_COLOR_SPACE i)
{
	if(i == OPJ_CLRSPC_SRGB) return "OPJ_CLRSPC_SRGB";
	if(i == OPJ_CLRSPC_GRAY) return "OPJ_CLRSPC_GRAY";
	if(i == OPJ_CLRSPC_SYCC) return "OPJ_CLRSPC_SYCC";
	if(i == OPJ_CLRSPC_UNKNOWN) return "OPJ_CLRSPC_UNKNOWN";
	return "CLRSPC_UNDEFINED";
}



void release(decode_info_t *decodeInfo)
{

	if(decodeInfo->codec)
	{
		opj_destroy_codec(decodeInfo->codec);
		decodeInfo->codec = NULL;
	}

	if(decodeInfo->stream)
	{
		opj_stream_destroy_v3(decodeInfo->stream);
		decodeInfo->stream = NULL;
	}

	if(decodeInfo->deleteImage && decodeInfo->image)
	{
		opj_image_destroy(decodeInfo->image);
		decodeInfo->image = NULL;

	}
}


OPJSupport::OPJSupport(){

}
OPJSupport::~OPJSupport(){

}


void* OPJSupport::decompressJPEG2K( void* jp2Data, long jp2DataSize, long *decompressedBufferSize, int *colorModel){

    return decompressJPEG2KWithBuffer(NULL, jp2Data, jp2DataSize, decompressedBufferSize, colorModel);


}
void* OPJSupport::decompressJPEG2KWithBuffer( void* inputBuffer, void* jp2Data, long jp2DataSize, long *decompressedBufferSize, int *colorModel){
    opj_dparameters_t parameters;
	OPJ_BOOL hasFile = OPJ_FALSE;

	opj_buffer_info_t buf_info;
	int i, decod_format;
	int width, height;
	OPJ_BOOL hasAlpha, fails = OPJ_FALSE;
	OPJ_CODEC_FORMAT codec_format;
	unsigned char rc, gc, bc, ac;

	decode_info_t decodeInfo;


	memset(&decodeInfo, 0, sizeof(decode_info_t));
	memset(&buf_info, 0, sizeof(opj_buffer_info_t));

	if (jp2Data != NULL)
	{
		buf_info.len = jp2DataSize;
		buf_info.buf =  (OPJ_BYTE*)jp2Data;
		buf_info.cur = buf_info.buf;
	}


	opj_set_default_decoder_parameters(&parameters);
	decod_format = buffer_format(&buf_info);

	if(decod_format == -1)
	{
		fprintf(stderr,"%s:%d: decode format missing\n",__FILE__,__LINE__);
		release(&decodeInfo);
		return 0;
	}

	/*-----------------------------------------------*/
	if(decod_format == J2K_CFMT)
		codec_format = OPJ_CODEC_J2K;
	else
		if(decod_format == JP2_CFMT)
			codec_format = OPJ_CODEC_JP2;
		else
			if(decod_format == JPT_CFMT)
				codec_format = OPJ_CODEC_JPT;
			else
			{
				/* clarified in infile_format() : */
				release(&decodeInfo);
				return 0;
			}
    parameters.decod_format = decod_format;
    while(1)
    {
	int tile_index=-1, user_changed_tile=0, user_changed_reduction=0;
	int max_tiles=0, max_reduction=0;
	fails = OPJ_TRUE;
	decodeInfo.stream =  opj_stream_create_buffer_stream(&buf_info, 1);


	if(decodeInfo.stream == NULL)
	{
	    fprintf(stderr,"%s:%d: NO decodeInfo.stream\n",__FILE__,__LINE__);
	    break;
	}
	decodeInfo.codec = opj_create_decompress(codec_format);
	if(decodeInfo.codec == NULL)
	{
	    fprintf(stderr,"%s:%d: NO coded\n",__FILE__,__LINE__);
	    break;
	}

      //  opj_set_info_handler(decodeInfo.codec, error_callback, this);
      //  opj_set_info_handler(decodeInfo.codec, warning_callback, this);
      //  opj_set_info_handler(decodeInfo.codec, info_callback, this);

	if( !opj_setup_decoder(decodeInfo.codec, &parameters))
	{
	    fprintf(stderr,"%s:%d:\n\topj_setup_decoder failed\n",__FILE__,__LINE__);
	    break;
	}

	if(user_changed_tile && user_changed_reduction)
	{
	    int reduction=0;
	    opj_set_decoded_resolution_factor(decodeInfo.codec, reduction);
	}

	if( !opj_read_header(decodeInfo.stream, decodeInfo.codec, &decodeInfo.image))
	{
	    fprintf(stderr,"%s:%d:\n\topj_read_header failed\n",__FILE__,__LINE__);
	    break;
	}

	if( !(user_changed_tile && user_changed_reduction)
	   || (max_tiles <= 0) || (max_reduction <= 0) )
	{
	    opj_codestream_info_v2_t *cstr;

	    cstr = opj_get_cstr_info(decodeInfo.codec);

	    max_reduction = cstr->m_default_tile_info.tccp_info->numresolutions;
	    max_tiles = cstr->tw * cstr->th;
	}

	if(tile_index < 0)
	{
	    unsigned int x0, y0, x1, y1;
	    int user_changed_area=0;

	    x0 = y0 = x1 = y1 = 0;


	    if(user_changed_area)
	    {

	    }

	    if( !opj_set_decode_area(decodeInfo.codec, decodeInfo.image, x0, y0, x1, y1))
	    {
		fprintf(stderr,"%s:%d:\n\topj_set_decode_area failed\n",__FILE__,__LINE__);
		break;
	    }
	    if( !opj_decode(decodeInfo.codec, decodeInfo.stream, decodeInfo.image))
	    {
		fprintf(stderr,"%s:%d:\n\topj_decode failed\n",__FILE__,__LINE__);
		break;
	    }
	}	/* if(tile_index < 0) */
	else
	{
	    if( !opj_get_decoded_tile(decodeInfo.codec, decodeInfo.stream, decodeInfo.image, tile_index))
	    {
		fprintf(stderr,"%s:%d:\n\topj_get_decoded_tile failed\n",__FILE__,__LINE__);
		break;
	    }
	}

	if( !opj_end_decompress(decodeInfo.codec, decodeInfo.stream))
	{
	    fprintf(stderr,"%s:%d:\n\topj_end_decompress failed\n",__FILE__,__LINE__);
	    break;
	}

	fails = OPJ_FALSE;
	break;

    }
    decodeInfo.deleteImage = fails;
    release(&decodeInfo);
    if(fails)
    {
	return 0;
    }
    decodeInfo.deleteImage = OPJ_TRUE;

    if(decodeInfo.image->color_space != OPJ_CLRSPC_SYCC
       && decodeInfo.image->numcomps == 3
       && decodeInfo.image->comps[0].dx == decodeInfo.image->comps[0].dy
       && decodeInfo.image->comps[1].dx != 1)
	decodeInfo.image->color_space = OPJ_CLRSPC_SYCC;
    else
	if(decodeInfo.image->numcomps <= 2)
	    decodeInfo.image->color_space = OPJ_CLRSPC_GRAY;

    if(decodeInfo.image->color_space == OPJ_CLRSPC_SYCC)
    {
	//disable for now
	//color_sycc_to_rgb(decodeInfo.image);
    }
    if(decodeInfo.image->icc_profile_buf)
    {
#if defined(HAVE_LIBLCMS1) || defined(HAVE_LIBLCMS2)
	color_apply_icc_profile(decodeInfo.image);
#endif

	free(decodeInfo.image->icc_profile_buf);
	decodeInfo.image->icc_profile_buf = NULL;
	decodeInfo.image->icc_profile_len = 0;
    }

    width = decodeInfo.image->comps[0].w;
    height = decodeInfo.image->comps[0].h;

    long depth = (decodeInfo.image->comps[0].prec + 7)/8;
    long decompressSize = width * height * decodeInfo.image->numcomps * depth;
    if (decompressedBufferSize)
	*decompressedBufferSize = decompressSize;;

    if (!inputBuffer ) {
	inputBuffer =  malloc(decompressSize);
    }

    if (colorModel)
	*colorModel = 0;

    if ((decodeInfo.image->numcomps >= 3
	 && decodeInfo.image->comps[0].dx == decodeInfo.image->comps[1].dx
	 && decodeInfo.image->comps[1].dx == decodeInfo.image->comps[2].dx
	 && decodeInfo.image->comps[0].dy == decodeInfo.image->comps[1].dy
	 && decodeInfo.image->comps[1].dy == decodeInfo.image->comps[2].dy
	 && decodeInfo.image->comps[0].prec == decodeInfo.image->comps[1].prec
	 && decodeInfo.image->comps[1].prec == decodeInfo.image->comps[2].prec
	 )/* RGB[A] */
	||
	(decodeInfo.image->numcomps == 2
	 && decodeInfo.image->comps[0].dx == decodeInfo.image->comps[1].dx
	 && decodeInfo.image->comps[0].dy == decodeInfo.image->comps[1].dy
	 && decodeInfo.image->comps[0].prec == decodeInfo.image->comps[1].prec
	 )
	) /* GA */
    {
	int  has_alpha4, has_alpha2, has_rgb;
	int *red, *green, *blue, *alpha;

	if (colorModel)
	    *colorModel = 1;

	alpha = NULL;

	has_rgb = (decodeInfo.image->numcomps == 3);
	has_alpha4 = (decodeInfo.image->numcomps == 4);
	has_alpha2 = (decodeInfo.image->numcomps == 2);
	hasAlpha = (has_alpha4 || has_alpha2);

	if(has_rgb)
	{
	    red = decodeInfo.image->comps[0].data;
	    green = decodeInfo.image->comps[1].data;
	    blue = decodeInfo.image->comps[2].data;

	    if(has_alpha4)
	    {
		alpha = decodeInfo.image->comps[3].data;
	    }

	}	/* if(has_rgb) */
	else
	{
	    red = green = blue = decodeInfo.image->comps[0].data;
	    if(has_alpha2)
	    {
		alpha = decodeInfo.image->comps[1].data;
	    }
	}	/* if(has_rgb) */


	ac = 255;/* 255: FULLY_OPAQUE; 0: FULLY_TRANSPARENT */


	int* ptrIBody = (int*)inputBuffer;
	for(i = 0; i < width*height; i++)
	{
	    rc = (unsigned char) *red++;
	    gc = (unsigned char)*green++;
	    bc = (unsigned char)*blue++;
	    if(hasAlpha)
	    {
		ac = (unsigned char)*alpha++;;
	    }
	    /*                         A        R          G       B
	     */
	    *ptrIBody++ = (int)((ac<<24) | (rc<<16) | (gc<<8) | bc);

	}	/* for(i) */
    }/* if (decodeInfo.image->numcomps >= 3  */
    else
	if(decodeInfo.image->numcomps == 1) /* Grey */
	{
	    /* 1 component 8 or 16 bpp decodeInfo.image
	     */
	    int *grey = decodeInfo.image->comps[0].data;
	    if(decodeInfo.image->comps[0].prec <= 8)
	    {

		char* ptrBBody = (char*)inputBuffer;
		for(i=0; i<width*height; i++)
		{
		    *ptrBBody++ = *grey++;
		}
		/* Replace image8 buffer:
		 */
	    }
	    else /* prec[9:16] */
	    {
		int *grey;
		int ushift = 0, dshift = 0, force16 = 0;

		grey = decodeInfo.image->comps[0].data;

		short* ptrSBody = (short*)inputBuffer;

		for(i=0; i<width*height; i++)
		{
		    //disable shift up for signed data: don't know why we are doing this
		    *ptrSBody++ = *grey++;
		}
		/* Replace image16 buffer:
		 */
	    }
	}
	else
	{
	    int *grey;

	    fprintf(stderr,"%s:%d:Can show only first component of decodeInfo.image\n"
		    "  components(%d) prec(%d) color_space[%d](%s)\n"
		    "  RECT(%d,%d,%d,%d)\n",__FILE__,__LINE__,decodeInfo.image->numcomps,
		    decodeInfo.image->comps[0].prec,
		    decodeInfo.image->color_space,clr_space(decodeInfo.image->color_space),
		    decodeInfo.image->x0,decodeInfo.image->y0,decodeInfo.image->x1,decodeInfo.image->y1 );

	    for(i = 0; i < decodeInfo.image->numcomps; ++i)
	    {
		fprintf(stderr,"[%d]dx(%d) dy(%d) w(%d) h(%d) signed(%u)\n",i,
			decodeInfo.image->comps[i].dx ,decodeInfo.image->comps[i].dy,
			decodeInfo.image->comps[i].w,decodeInfo.image->comps[i].h,
			decodeInfo.image->comps[i].sgnd);
	    }

	    /* 1 component 8 or 16 bpp decodeInfo.image
	     */
	    grey = decodeInfo.image->comps[0].data;
	    if(decodeInfo.image->comps[0].prec <= 8)
	    {

		char* ptrBBody = (char*)inputBuffer;
		for(i=0; i<width*height; i++)
		{
		    *ptrBBody++ = *grey++;
		}
		/* Replace image8 buffer:
		 */
	    }
	    else /* prec[9:16] */
	    {
		int *grey;
		int ushift = 0, dshift = 0, force16 = 0;

		grey = decodeInfo.image->comps[0].data;

		short* ptrSBody = (short*)inputBuffer;

		for(i=0; i<width*height; i++)
		{
		    *ptrSBody++ = *grey++;
		}
		/* Replace image16 buffer:
		 */
	    }
	}
    release(&decodeInfo);
     return inputBuffer;
}
void* OPJSupport::compressJPEG2K( void *data, int samplesPerPixel, int rows, int columns, int precision, bool sign, int rate, long *compressedDataSize){
    return 0;
}
