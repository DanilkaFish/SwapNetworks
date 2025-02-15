import express, { NextFunction } from "express";
import { Request, Response } from "express";
import { jwtVerify } from "../util/jwt.utils";
import { decode } from "punycode";
import { throwError } from "../util/utils";
// export const authRouter = express.Router();


// authRouter.post("/login", userContoller.login);
// authRouter.post("/refresh_token", userContoller.refreshToken);
// authRouter.post("/logout", userContoller.logout);


export const auth = async (req: Request, res: Response, next: NextFunction): Promise<void> => {
  try {
    const token = req.header("Authorization")?.replace("Bearer ", "");
    const { email } = req.body;
    if (!token) {
      res.status(401).json({ message: 'Unauthorized' });
    } else {
      const decoded = jwtVerify(token);
      if (decoded){
        // req.userid = 
        next()
      } else{
        throwError(
          "Unauthorized",
          401,
          "Unauthorized"
        );
      }
    }
    } catch (err) {
      res.status(401).json({ message: 'Invalid Token' });
  }
};
