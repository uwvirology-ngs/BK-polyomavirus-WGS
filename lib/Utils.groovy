class Utils {
    public static String getAnnotation(String s) {
        String ref = s.substring(s.lastIndexOf('_') + 1)
        return ref
    }

    public static String getGenomicRegion(String s) {
        HashMap map = [
            "PZ129802.1": "PZ129802.1:1470-2558",
        ]
        return map[s]
    }
}