import { Application, Request, Response } from "express";
import { NextFunction } from "express";
import { errorResponce } from "../util/utils";
import { authRouter } from "./auth.routes";

export const setupRoutes = (app: Application): void => {
    app.get("/", (req: Request, res: Response) => {
        res.send("Hello World! Alles!");
    });
    app.use("/authorization", authRouter);
};
