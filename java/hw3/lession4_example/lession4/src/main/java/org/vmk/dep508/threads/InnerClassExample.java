package org.vmk.dep508.threads;

/**
 * Created by Asus on 10/1/2018.
 */
public class InnerClassExample {
    public static void main(String[] args) {

        InnerClassExample example = new InnerClassExample();
        example.innerClassInClass();
        example.innerClsssInMethod();

        StaticLogger.log("Hello from static class");

        StaticLogger logger = new StaticLogger();
        logger.nonStaticLog("Hello from static class");

    }

    String var = "some class variable";
    class ClassLogger {
        String name;

        public ClassLogger(String name) {
            this.name = name;
            System.out.println(var);
        }

        public void log(String message) {

            System.out.println("class "+name + " - "+message);
        }
    }

    public void innerClassInClass() {
        ClassLogger logger = new ClassLogger("innerClassInClass");
        logger.log("Hello from innerClassInClass");
    }

    static class StaticLogger {
        static void log(String msg) {
            System.out.println("static inner class - "+msg);
        }
        void nonStaticLog(String msg) {
            System.out.println("static class, non static method - "+msg);
        }
    }

    public void innerClsssInMethod() {

        class MethodLogger {
            String name;

            public MethodLogger(String name) {
                this.name = name;
            }

            public void log(String message) {
                System.out.println(name+" - "+message);
            }
        }

        MethodLogger logger = new MethodLogger("innerClsssInMethodLogger");

        logger.log("Hello from inner class");
    }
}
