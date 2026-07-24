class Utils {
    public static String getAnnotation(String s) {
        String ref = s.substring(s.indexOf('_') + 1)
        return ref
    }

    public static String getGenomicRegion(String s) {
        HashMap map = [
            "PV856408.1": "PV856408.1:1177-2265",
            "PV856409.1": "PV856409.1:1177-2265",
            "PV856429.1": "PV856429.1:1167-2255",
            "PV856444.1": "PV856444.1:1177-2265",
            "PV856452.1": "PV856452.1:1177-2265",
            "PX257737.1": "PX257737.1:1193-2281",
            "PZ129835.1": "PZ129835.1:1470-2558",
            "NC_001538.1": "NC_001538.1:1564-2652"
        ]
        return map[s]
    }
}