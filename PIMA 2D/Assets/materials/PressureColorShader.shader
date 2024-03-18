Shader "Custom/PressureColorShader"
{
    Properties
    {
        _Pressure ("Pressure", Range(0, 1)) = 0.5
    }
    SubShader
    {
        Tags { "RenderType"="Opaque" }
        LOD 100

        Pass
        {
            CGPROGRAM
            #pragma vertex vert
            #pragma fragment frag

            #include "UnityCG.cginc"

            struct appdata
            {
                float4 vertex : POSITION;
                float2 uv : TEXCOORD0;
            };

            struct v2f
            {
                float2 uv : TEXCOORD0;
                UNITY_FOG_COORDS(1)
                float4 vertex : SV_POSITION;
            };

            float _Pressure;

            v2f vert (appdata v)
            {
                v2f o;
                o.vertex = UnityObjectToClipPos(v.vertex);
                o.uv = v.uv;
                return o;
            }

            fixed4 frag (v2f i) : SV_Target
            {
                // Define colors for high and low pressure
                fixed4 redColor = fixed4(1, 0, 0, 1);
                fixed4 blueColor = fixed4(0, 0, 1, 1);

                // Interpolate between red and blue based on pressure value
                fixed4 finalColor = lerp(blueColor, redColor, _Pressure);

                return finalColor;
            }
            ENDCG
        }
    }
}